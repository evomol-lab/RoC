"""The Room of Conformations - web interface (Flask).

Every calculation runs roc.py in its own subprocess and job directory, so long
ensembles can be followed (progress + log) and cancelled from the browser.

    python RoCGUI.py                 # http://127.0.0.1:5050
    ROC_HOST=0.0.0.0 python RoCGUI.py   # listen on all interfaces (containers)
"""
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import uuid

from flask import Flask, abort, jsonify, render_template, request, send_file, send_from_directory
from werkzeug.utils import secure_filename

FROZEN = getattr(sys, 'frozen', False)
APP_ROOT = sys._MEIPASS if FROZEN else os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, template_folder=os.path.join(APP_ROOT, 'templates'),
            static_folder=os.path.join(APP_ROOT, 'static'))
app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024

JOBS_DIR = os.environ.get('ROC_JOBS_DIR', os.path.abspath('roc_jobs'))
os.makedirs(JOBS_DIR, exist_ok=True)
RETENTION_DAYS = float(os.environ.get('ROC_JOB_RETENTION_DAYS', '0'))   # 0 keeps every job

STRUCTURE_EXT = ('.pdb', '.ent', '.cif', '.mmcif', '.pdb.gz', '.ent.gz', '.cif.gz', '.mmcif.gz')
JOB_ID = re.compile(r'^[0-9a-f]{12}$')
PROCS = {}   # job id -> subprocess.Popen


def clustenm_available():
    return all(importlib.util.find_spec(m) is not None for m in ('openmm', 'pdbfixer'))


def job_dir(job_id):
    if not JOB_ID.match(job_id or ''):
        abort(404)
    path = os.path.join(JOBS_DIR, job_id)
    if not os.path.isdir(path):
        abort(404)
    return path


def purge_old_jobs():
    if RETENTION_DAYS <= 0:
        return
    limit = time.time() - RETENTION_DAYS * 86400
    for name in os.listdir(JOBS_DIR):
        path = os.path.join(JOBS_DIR, name)
        if JOB_ID.match(name) and os.path.isdir(path) and os.path.getmtime(path) < limit:
            shutil.rmtree(path, ignore_errors=True)


def fetch_structure(db, ident, dest_dir):
    """Download a structure from RCSB PDB or AlphaFold DB; returns the saved path."""
    ident = ident.strip()
    if db == 'pdb':
        if not re.match(r'^[0-9][A-Za-z0-9]{3}$', ident):
            raise ValueError('A PDB ID has four characters, e.g. 4AKE.')
        url = f'https://files.rcsb.org/download/{ident.upper()}.cif'
        name = f'{ident.upper()}.cif'
    elif db == 'alphafold':
        if not re.match(r'^[A-Za-z0-9]{6,10}$', ident):
            raise ValueError('Use a UniProt accession, e.g. P61851.')
        acc = ident.upper()
        url = None
        try:
            with urllib.request.urlopen(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}',
                                        timeout=30) as r:
                entries = json.loads(r.read().decode())
            if entries:
                url = entries[0].get('cifUrl') or entries[0].get('pdbUrl')
        except Exception:
            pass
        url = url or f'https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif'
        name = os.path.basename(url.split('?')[0]) or f'AF-{acc}.cif'
    else:
        raise ValueError('Unknown database.')
    path = os.path.join(dest_dir, secure_filename(name))
    try:
        with urllib.request.urlopen(url, timeout=60) as r, open(path, 'wb') as fh:
            shutil.copyfileobj(r, fh)
    except Exception as exc:
        raise ValueError(f'Could not download {ident} ({exc}).')
    return path


def save_upload(storage, dest_dir):
    name = secure_filename(storage.filename or '')
    if not name.lower().endswith(STRUCTURE_EXT):
        raise ValueError('Upload a .pdb, .ent or .cif file (optionally gzipped).')
    path = os.path.join(dest_dir, name)
    storage.save(path)
    return path


def number(form, key, cast, default, lo, hi):
    raw = form.get(key, '')
    if raw in ('', None):
        return default
    try:
        value = cast(raw)
    except ValueError:
        raise ValueError(f'Invalid value for {key}: {raw}')
    if not lo <= value <= hi:
        raise ValueError(f'{key} must be between {lo} and {hi}.')
    return value


def build_arguments(form, action):
    """Validated roc.py command-line options from the form."""
    choice = lambda key, allowed, default: form.get(key) if form.get(key) in allowed else default
    method = choice('method', ('anm', 'rtb', 'clustenm'), 'anm')
    if method == 'clustenm' and not clustenm_available():
        raise ValueError('ClustENM needs OpenMM and PDBFixer, which are not installed.')
    args = ['--method', method,
            '--modes', str(number(form, 'modes', int, 15, 1, 100)),
            '--confs', str(number(form, 'confs', int, 50, 1, 5000)),
            '--rmsd', str(number(form, 'rmsd', float, 1.0, 0.05, 20)),
            '--sampling', choice('sampling', ('random', 'traverse'), 'random'),
            '--gamma', str(number(form, 'gamma', float, 1.0, 0.001, 1000)),
            '--blocks', choice('blocks', ('residue', 'sse'), 'residue'),
            '--springs', choice('springs', ('uniform', 'structure'), 'uniform'),
            '--ligands', choice('ligands', ('follow', 'network', 'remove'), 'follow'),
            '--clusters', str(number(form, 'clusters', int, 5, 1, 50)),
            '--trim-plddt', str(number(form, 'trim_plddt', float, 0, 0, 100)),
            '--clustenm-gens', str(number(form, 'clustenm_gens', int, 1, 1, 5))]
    cutoff = number(form, 'cutoff', float, None, 3, 50)
    if cutoff is not None:
        args += ['--cutoff', str(cutoff)]
    seed = number(form, 'seed', int, None, 0, 2 ** 31 - 1)
    if seed is not None:
        args += ['--seed', str(seed)]
    chains = form.get('chains', '').strip()
    if chains:
        if not re.match(r'^[A-Za-z0-9]{1,4}(\s*,\s*[A-Za-z0-9]{1,4})*$', chains):
            raise ValueError('Chains must be IDs separated by commas, e.g. A,B')
        args += ['--chains', chains.replace(' ', '')]
    if form.get('split_models') == 'true':
        args.append('--split-models')
    if form.get('clustenm_md') == 'true':
        args.append('--clustenm-md')
    if action == 'analyze':
        args.append('--analyze-only')
    return args


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/documentation')
def documentation():
    path = os.path.join(APP_ROOT, 'DOCUMENTATION.md')
    if os.path.exists(path):
        return send_file(path, mimetype='text/markdown')
    return 'Documentation not found', 404


@app.route('/api/capabilities')
def capabilities():
    import roc
    return jsonify({'version': roc.__version__, 'clustenm': clustenm_available(),
                    'analysis_modes': roc.N_ANALYSIS_MODES,
                    'default_cutoff': roc.DEFAULT_CUTOFF})


@app.route('/api/jobs', methods=['POST'])
def submit():
    form = request.form
    action = 'analyze' if form.get('action') == 'analyze' else 'ensemble'
    job_id = uuid.uuid4().hex[:12]
    path = os.path.join(JOBS_DIR, job_id)
    inputs = os.path.join(path, 'input')          # inputs are kept apart from the results
    target_dir = os.path.join(inputs, 'target')
    os.makedirs(target_dir)
    target = None
    try:
        args = build_arguments(form, action)
        mode = form.get('input_mode', 'upload')
        if mode == 'reuse':
            src = os.path.join(JOBS_DIR, form.get('reuse_job', ''))
            if not JOB_ID.match(form.get('reuse_job', '')) or not os.path.isdir(src):
                raise ValueError('The previous job was not found; upload the structure again.')
            with open(os.path.join(src, 'job.json')) as fh:
                old = json.load(fh)
            structure = shutil.copy(os.path.join(src, 'input', old['input']), inputs)
            if old.get('target') and form.get('keep_target') != 'false':
                target = shutil.copy(os.path.join(src, 'input', 'target', old['target']),
                                     target_dir)
        elif mode == 'fetch':
            structure = fetch_structure(form.get('fetch_db'), form.get('fetch_id', ''), inputs)
        else:
            if 'file' not in request.files or not request.files['file'].filename:
                raise ValueError('No structure file was uploaded.')
            structure = save_upload(request.files['file'], inputs)
        if 'target' in request.files and request.files['target'].filename:
            target = save_upload(request.files['target'], target_dir)
    except ValueError as exc:
        shutil.rmtree(path, ignore_errors=True)
        return jsonify({'error': str(exc)}), 400
    if target:
        args += ['--target', target]
    name = os.path.splitext(os.path.basename(structure))[0]
    name = re.sub(r'\.(pdb|cif|ent|mmcif)$', '', name, flags=re.I)
    args += ['--name', name]
    with open(os.path.join(path, 'job.json'), 'w') as fh:
        json.dump({'id': job_id, 'created': time.time(), 'action': action,
                   'input': os.path.basename(structure),
                   'target': os.path.basename(target) if target else None, 'args': args}, fh)

    if FROZEN:
        cmd = [sys.executable, '--run-roc', structure, '--outdir', path] + args
    else:
        cmd = [sys.executable, os.path.join(APP_ROOT, 'roc.py'), structure, '--outdir', path] + args
    log = open(os.path.join(path, 'run.log'), 'w')
    env = dict(os.environ, PYTHONUNBUFFERED='1', MPLBACKEND='Agg')
    PROCS[job_id] = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=path, env=env)
    log.close()
    return jsonify({'job': job_id})


@app.route('/api/jobs/<job_id>')
def job_status(job_id):
    path = job_dir(job_id)
    proc = PROCS.get(job_id)
    status_file = os.path.join(path, 'status.json')
    state, error = 'running', None
    if os.path.exists(status_file):
        with open(status_file) as fh:
            st = json.load(fh)
        state, error = st['state'], st.get('error')
    elif proc is None or proc.poll() is not None:
        if os.path.exists(os.path.join(path, 'cancelled')):
            state = 'cancelled'
        else:
            state, error = 'failed', 'The calculation stopped unexpectedly (see the log).'
    lines, prog, message = [], 0.0, ''
    log_path = os.path.join(path, 'run.log')
    if os.path.exists(log_path):
        with open(log_path, errors='replace') as fh:
            for line in fh:
                if line.startswith('@@progress'):
                    parts = line.split(maxsplit=2)
                    prog = float(parts[1])
                    message = parts[2].strip() if len(parts) > 2 else ''
                elif 'Warning' not in line and not line.startswith('  warnings.warn'):
                    lines.append(line.rstrip('\n'))
    out = {'state': state, 'error': error, 'progress': prog, 'message': message,
           'log': '\n'.join(lines[-400:])}
    if state == 'done':
        with open(os.path.join(path, 'summary.json')) as fh:
            out['summary'] = json.load(fh)
    return jsonify(out)


@app.route('/api/jobs/<job_id>/cancel', methods=['POST'])
def cancel(job_id):
    path = job_dir(job_id)
    proc = PROCS.get(job_id)
    if proc and proc.poll() is None:
        open(os.path.join(path, 'cancelled'), 'w').close()
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
    return jsonify({'state': 'cancelled'})


@app.route('/api/jobs/<job_id>/files/<path:filename>')
def job_file(job_id, filename):
    path = job_dir(job_id)
    return send_from_directory(path, filename, as_attachment=request.args.get('download') == '1')


@app.route('/api/jobs/<job_id>/mode/<int:mode>')
def mode_animation(job_id, mode):
    path = job_dir(job_id)
    target = os.path.join(path, f'mode_{mode:03d}.pdb')
    if not os.path.exists(target):
        import roc
        try:
            roc.mode_trajectory(path, mode)
        except roc.RocError as exc:
            return jsonify({'error': str(exc)}), 400
    return send_file(target, mimetype='chemical/x-pdb')


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--run-roc':
        import roc
        sys.exit(roc.main(sys.argv[2:]))
    purge_old_jobs()
    host = os.environ.get('ROC_HOST', '127.0.0.1')
    port = int(os.environ.get('ROC_PORT', '5050'))
    url = f'http://{"127.0.0.1" if host == "0.0.0.0" else host}:{port}/'
    print(f'The Room of Conformations: open {url} in your browser (Ctrl+C to stop).')
    print(f'Jobs are stored in {JOBS_DIR}')
    if FROZEN or '--browser' in sys.argv:
        import threading
        import webbrowser
        threading.Timer(1.25, lambda: webbrowser.open_new(url)).start()
    # Debug mode stays off: Werkzeug's interactive debugger would let anyone who
    # reaches the port run code.
    app.run(host=host, port=port, debug=False, threaded=True)


if __name__ == '__main__':
    main()
