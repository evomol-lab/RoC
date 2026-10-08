'use strict';

document.addEventListener('DOMContentLoaded', () => {
    const $ = (id) => document.getElementById(id);
    const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({
        '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    })[c]);
    const fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v)) ? '–' : Number(v).toFixed(d);
    const pct = (v, d = 1) => (v === null || v === undefined) ? '–' : (100 * v).toFixed(d) + '%';

    const state = {
        caps: null,
        file: null,
        targetFile: null,
        job: null,           // job being run
        resultJob: null,     // job whose results are shown
        inputJob: null,      // last job whose input can be reused
        inputDirty: true,
        targetDirty: false,
        summary: null,
        K: 15,
        modeViewer: null,
        modePlaying: true,
        ensViewer: null,
        ensTimer: null,
        ensFrame: 0,
        ensFrames: 0,
        pollTimer: null,
    };

    // ------------------------------------------------------------------ theme
    const themeIconLight = $('theme-icon-light');
    const themeIconDark = $('theme-icon-dark');

    function storedTheme() {
        try { return localStorage.getItem('roc-theme'); } catch (e) { return null; }
    }

    function applyTheme(light, persist) {
        document.body.classList.toggle('light-mode', light);
        themeIconDark.style.display = light ? 'none' : 'block';
        themeIconLight.style.display = light ? 'block' : 'none';
        if (persist) {
            try { localStorage.setItem('roc-theme', light ? 'light' : 'dark'); } catch (e) { /* private mode */ }
        }
        if (state.summary) renderCharts();
        for (const v of [state.modeViewer, state.ensViewer]) if (v) {
            v.setBackgroundColor(viewerBackground(), 1);
            recolorViewer(v);
            v.render();
        }
    }

    function viewerBackground() {
        return document.body.classList.contains('light-mode') ? '#ffffff' : '#14161b';
    }

    const saved = storedTheme();
    applyTheme(saved ? saved === 'light'
        : window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches, false);
    $('theme-toggle').addEventListener('click', () =>
        applyTheme(!document.body.classList.contains('light-mode'), true));

    // ------------------------------------------------------------------ documentation tab
    let docsLoaded = false;
    function showView(view) {
        const docs = view === 'documentation';
        $('nav-calculator').classList.toggle('active', !docs);
        $('nav-documentation').classList.toggle('active', docs);
        $('calculator-view').style.display = docs ? 'none' : 'block';
        $('documentation-view').style.display = docs ? 'block' : 'none';
        document.querySelector('.container').classList.toggle('expanded', docs || !!state.summary);
        if (docs && !docsLoaded) {
            fetch('/documentation').then((r) => {
                if (!r.ok) throw new Error('Documentation file not found on the server.');
                return r.text();
            }).then((md) => {
                $('doc-content').innerHTML = marked.parse(md);
                docsLoaded = true;
            }).catch((e) => {
                $('doc-content').innerHTML = `<div class="doc-loading">${esc(e.message)}</div>`;
            });
        }
        if (!docs && state.summary) setTimeout(resizeCharts, 50);
    }
    $('nav-calculator').addEventListener('click', () => showView('calculator'));
    $('nav-documentation').addEventListener('click', () => showView('documentation'));

    // ------------------------------------------------------------------ capabilities
    fetch('/api/capabilities').then((r) => r.json()).then((caps) => {
        state.caps = caps;
        const opt = document.querySelector('#method option[value="clustenm"]');
        if (!caps.clustenm) {
            opt.disabled = true;
            opt.textContent += ' (needs OpenMM + PDBFixer)';
        }
        updateMethod();
    }).catch(() => updateMethod());

    // ------------------------------------------------------------------ input
    const dropZone = $('drop-zone');
    const fileInput = $('file-input');
    const STRUCT_RE = /\.(pdb|ent|cif|mmcif)(\.gz)?$/i;

    function inputMode() {
        return document.querySelector('input[name="input_mode"]:checked').value;
    }

    document.querySelectorAll('input[name="input_mode"]').forEach((el) => el.addEventListener('change', () => {
        const fetchMode = inputMode() === 'fetch';
        dropZone.style.display = fetchMode ? 'none' : 'block';
        $('fetch-zone').style.display = fetchMode ? 'block' : 'none';
        state.inputDirty = true;
        updateButtons();
    }));

    function setFile(file) {
        if (!file) return;
        if (!STRUCT_RE.test(file.name)) {
            $('file-name').textContent = 'Unsupported file: use .pdb, .ent or .cif (optionally .gz)';
            state.file = null;
        } else {
            state.file = file;
            $('file-name').textContent = file.name;
            $('error-message').style.display = 'none';
        }
        state.inputDirty = true;
        updateButtons();
    }

    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach((ev) => dropZone.addEventListener(ev, (e) => {
        e.preventDefault();
        e.stopPropagation();
    }));
    ['dragenter', 'dragover'].forEach((ev) => dropZone.addEventListener(ev, () => dropZone.classList.add('dragover')));
    ['dragleave', 'drop'].forEach((ev) => dropZone.addEventListener(ev, () => dropZone.classList.remove('dragover')));
    dropZone.addEventListener('drop', (e) => setFile(e.dataTransfer.files[0]));
    dropZone.addEventListener('click', () => fileInput.click());
    dropZone.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); fileInput.click(); }
    });
    fileInput.addEventListener('change', () => setFile(fileInput.files[0]));
    $('fetch_id').addEventListener('input', () => { state.inputDirty = true; updateButtons(); });
    $('fetch_db').addEventListener('change', () => {
        $('fetch_id').placeholder = $('fetch_db').value === 'pdb' ? 'e.g. 4AKE' : 'e.g. P61851';
        state.inputDirty = true;
        updateButtons();
    });
    $('target-input').addEventListener('change', () => {
        state.targetFile = $('target-input').files[0] || null;
        state.targetDirty = true;
    });

    function hasInput() {
        if (state.inputJob && !state.inputDirty) return true;
        return inputMode() === 'fetch' ? $('fetch_id').value.trim() !== '' : !!state.file;
    }

    function running() {
        return !!state.pollTimer;
    }

    function updateButtons() {
        const ok = hasInput() && !running();
        $('analyze-btn').disabled = !ok;
        $('generate-btn').disabled = !ok;
    }

    // ------------------------------------------------------------------ parameters
    const METHOD_DESC = {
        anm: 'One node per C-alpha atom (Bahar lab ANM, as in ProDy). Fast; every residue is rebuilt as a rigid body oriented by its neighbours.',
        rtb: 'All heavy atoms are nodes; residues (or helices and strands) move as rigid blocks (Tama et al. 2000). Respects atomic packing: far fewer steric overlaps in compact proteins. Recommended for docking ensembles.',
        clustenm: 'ProDy ClustENM (Kurkcuoglu, Bahar & Doruker 2016): ANM sampling, clustering and OpenMM energy minimisation with hydrogens added. Best stereochemistry, but minutes to hours and ligands are removed.',
    };

    function updateMethod() {
        const m = $('method').value;
        $('method-desc').textContent = METHOD_DESC[m];
        $('springs-field').style.display = m === 'anm' ? 'flex' : 'none';
        $('blocks-field').style.display = m === 'rtb' ? 'flex' : 'none';
        $('clustenm-gens-field').style.display = m === 'clustenm' ? 'flex' : 'none';
        $('clustenm-md-field').style.display = m === 'clustenm' ? 'flex' : 'none';
        $('cutoff').placeholder = m === 'rtb' ? '8' : '15';
        $('sampling').disabled = m === 'clustenm';
        if (m === 'clustenm') $('sampling').value = 'random';
        document.querySelectorAll('input[name="ligands"]').forEach((el) => { el.disabled = m === 'clustenm'; });
        updateSampling();
    }
    $('method').addEventListener('change', updateMethod);

    function updateSampling() {
        const traverse = $('sampling').value === 'traverse';
        $('confs-label').textContent = traverse ? 'Frames per mode' : ($('method').value === 'clustenm'
            ? 'Conformers' : 'Conformations');
        $('rmsd-label').textContent = traverse ? 'C-alpha RMSD at both ends of each mode (Å)'
            : 'Average C-alpha RMSD to the input (Å)';
        if (traverse && Number($('confs').value) > 41) $('confs').value = 11;
    }
    $('sampling').addEventListener('change', updateSampling);

    // number of modes: slider + box
    const modesRange = $('modes');
    const modesNum = $('modes-num');
    let rafPending = false;

    function setK(k, fromUser = true) {
        const max = Number(modesRange.max);
        k = Math.max(1, Math.min(max, Math.round(Number(k) || 1)));
        state.K = k;
        modesRange.value = k;
        if (Number(modesNum.value) !== k) modesNum.value = k;
        if (!rafPending) {
            rafPending = true;
            requestAnimationFrame(() => {
                rafPending = false;
                updateModesLive();
                if (state.summary && fromUser) updateKDependent();
            });
        }
    }
    modesRange.addEventListener('input', () => setK(modesRange.value));
    modesNum.addEventListener('change', () => setK(modesNum.value));

    // ------------------------------------------------------------------ mode statistics
    function modeStats(K) {
        const s = state.summary;
        if (!s) return null;
        const m = s.modes;
        const n = m.eigvals.length;
        K = Math.min(K, n);
        let share = 0, coll = 0;
        const localized = [];
        for (let k = 0; k < K; k++) {
            share += m.variance_fraction[k];
            coll += m.collectivity[k];
            if (m.collectivity[k] < 0.15) localized.push(k + 1);
        }
        const out = { K, n, share, coll: coll / K, localized, exact: m.exact };
        if (s.target) {
            out.cumOverlap = s.target.cumulative_overlap[K - 1];
            out.residual = s.target.residual_rmsd[K - 1];
            out.targetRmsd = s.target.rmsd;
            const need = (thr) => {
                const i = s.target.cumulative_overlap.findIndex((v) => v >= thr);
                return i < 0 ? null : i + 1;
            };
            out.need70 = need(0.7);
            out.need90 = need(0.9);
        }
        // modes needed to reach shares of the fluctuation
        let acc = 0;
        out.need = {};
        for (let k = 0; k < n; k++) {
            acc += m.variance_fraction[k];
            for (const t of [0.25, 0.5, 0.75]) if (!out.need[t] && acc >= t) out.need[t] = k + 1;
        }
        return out;
    }

    function updateModesLive() {
        const st = modeStats(state.K);
        const box = $('modes-live');
        if (!st) {
            box.innerHTML = `The <strong>${state.K}</strong> slowest (lowest-frequency) modes will be combined. ` +
                'Run <em>Analyze modes</em> to see how much of the network\'s motion they capture and which motions they are.';
            return;
        }
        let html = `Modes <strong>1–${st.K}</strong> carry <strong>${pct(st.share)}</strong> of the network's ` +
            `fluctuation${st.exact ? '' : ' (of the computed modes)'}; mean collectivity <strong>${fmt(st.coll)}</strong>.`;
        if (st.localized.length) {
            html += ` ${st.localized.length} of them ${st.localized.length === 1 ? 'is' : 'are'} localised ` +
                `(mode${st.localized.length === 1 ? '' : 's'} ${st.localized.slice(0, 8).join(', ')}${st.localized.length > 8 ? '…' : ''}).`;
        }
        if (st.cumOverlap !== undefined) {
            html += ` They describe <strong>${fmt(st.cumOverlap)}</strong> of the change to the second conformation ` +
                `(cumulative overlap); RMSD ${fmt(st.targetRmsd)} → ${fmt(st.residual)} Å.`;
        }
        box.innerHTML = html;
    }

    // ------------------------------------------------------------------ submission & polling
    function collectParams(fd) {
        fd.append('method', $('method').value);
        fd.append('modes', state.K);
        fd.append('confs', $('confs').value);
        fd.append('rmsd', $('rmsd').value);
        fd.append('sampling', $('sampling').value);
        fd.append('cutoff', $('cutoff').value);
        fd.append('gamma', $('gamma').value);
        fd.append('blocks', $('blocks').value);
        fd.append('springs', $('springs').value);
        fd.append('ligands', document.querySelector('input[name="ligands"]:checked').value);
        fd.append('chains', $('chains').value);
        fd.append('trim_plddt', $('trim_plddt').value);
        fd.append('seed', $('seed').value);
        fd.append('clusters', $('clusters').value);
        fd.append('split_models', $('split_models').checked);
        fd.append('clustenm_gens', $('clustenm_gens').value);
        fd.append('clustenm_md', $('clustenm_md').checked);
    }

    function showError(msg) {
        const box = $('error-message');
        box.textContent = msg;
        box.style.display = msg ? 'block' : 'none';
    }

    async function submit(action) {
        showError('');
        const fd = new FormData();
        fd.append('action', action);
        if (state.inputJob && !state.inputDirty) {
            fd.append('input_mode', 'reuse');
            fd.append('reuse_job', state.inputJob);
            fd.append('keep_target', state.targetDirty ? 'false' : 'true');
            if (state.targetDirty && state.targetFile) fd.append('target', state.targetFile);
        } else {
            const mode = inputMode();
            fd.append('input_mode', mode);
            if (mode === 'fetch') {
                fd.append('fetch_db', $('fetch_db').value);
                fd.append('fetch_id', $('fetch_id').value.trim());
            } else {
                fd.append('file', state.file);
            }
            if (state.targetFile) fd.append('target', state.targetFile);
        }
        collectParams(fd);
        $('analyze-btn').disabled = $('generate-btn').disabled = true;
        $('generate-btn').classList.add('loading');
        try {
            const r = await fetch('/api/jobs', { method: 'POST', body: fd });
            const d = await r.json();
            if (!r.ok) throw new Error(d.error || 'The job could not be started.');
            state.job = d.job;
            state.inputJob = d.job;
            state.inputDirty = false;
            state.targetDirty = false;
            state.jobAction = action;
            $('job-panel').style.display = 'block';
            $('progress-bar').style.width = '0%';
            $('progress-text').textContent = action === 'analyze' ? 'Analysing normal modes…' : 'Generating the ensemble…';
            $('log-box').textContent = '';
            state.pollTimer = setTimeout(poll, 400);
            history.replaceState(null, '', `?job=${d.job}`);
        } catch (e) {
            showError(e.message);
            $('generate-btn').classList.remove('loading');
        }
        updateButtons();
    }

    $('analyze-btn').addEventListener('click', () => submit('analyze'));
    $('generate-btn').addEventListener('click', () => submit('ensemble'));

    async function poll() {
        let d;
        try {
            const r = await fetch(`/api/jobs/${state.job}`);
            d = await r.json();
        } catch (e) {
            state.pollTimer = setTimeout(poll, 2000);
            return;
        }
        $('progress-bar').style.width = `${Math.round(100 * d.progress)}%`;
        if (d.message) $('progress-text').textContent = d.message + '…';
        const logBox = $('log-box');
        const atBottom = logBox.scrollTop + logBox.clientHeight >= logBox.scrollHeight - 20;
        logBox.textContent = d.log;
        if (atBottom) logBox.scrollTop = logBox.scrollHeight;
        if (d.state === 'running') {
            state.pollTimer = setTimeout(poll, 1000);
            return;
        }
        state.pollTimer = null;
        $('generate-btn').classList.remove('loading');
        updateButtons();
        if (d.state === 'done') {
            $('progress-text').textContent = 'Done.';
            setTimeout(() => { if (!running()) $('job-panel').style.display = 'none'; }, 1200);
            state.resultJob = state.job;
            showResults(d.summary);
        } else if (d.state === 'cancelled') {
            $('progress-text').textContent = 'Cancelled.';
        } else {
            $('progress-text').textContent = 'Failed.';
            logBox.parentElement.open = true;
            showError(d.error || 'The calculation failed; see the log.');
        }
    }

    $('cancel-btn').addEventListener('click', async () => {
        if (!state.job || !running()) return;
        await fetch(`/api/jobs/${state.job}/cancel`, { method: 'POST' });
    });

    // ------------------------------------------------------------------ results
    const fileUrl = (name, download) =>
        `/api/jobs/${state.resultJob}/files/${encodeURIComponent(name)}${download ? '?download=1' : ''}`;
    const METHOD_NAME = { anm: 'ANM', rtb: 'RTB', clustenm: 'ClustENM' };

    function showResults(summary) {
        state.summary = summary;
        const n = summary.modes.eigvals.length;
        modesRange.max = Math.min(100, n);
        modesNum.max = modesRange.max;
        setK(summary.parameters.modes, false);
        $('results').style.display = 'block';
        document.querySelector('.container').classList.add('expanded');
        const p = summary.parameters;
        $('results-title').textContent = `${summary.name} · ${METHOD_NAME[p.method]}` +
            (summary.ensemble ? ` · ${summary.ensemble.n_models} models` : ' · normal modes');
        $('zip-btn').href = fileUrl(summary.files.zip, true);
        const hasEns = !!summary.ensemble;
        $('tab-ensemble-btn').disabled = !hasEns;
        switchTab(hasEns ? 'ensemble' : 'modes');
        renderModesStatic();
        renderCharts();
        updateKDependent();
        renderFiles();
        fillModeSelect();
        showMode(1);
        if (hasEns) {
            renderEnsembleStatic();
            loadEnsembleViewer();
        } else if (state.ensViewer) {
            stopEnsemble();
            state.ensViewer.clear();
            state.ensViewer.render();
        }
        updateModesLive();
        $('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    function switchTab(name) {
        document.querySelectorAll('.tab-btn').forEach((b) => b.classList.toggle('active', b.dataset.tab === name));
        document.querySelectorAll('.tab-pane').forEach((p) => p.classList.toggle('active', p.id === `pane-${name}`));
        setTimeout(() => {
            resizeCharts();
            for (const v of [state.modeViewer, state.ensViewer]) if (v) { v.resize(); v.render(); }
        }, 30);
    }
    document.querySelectorAll('.tab-btn').forEach((b) => b.addEventListener('click', () => {
        if (!b.disabled) switchTab(b.dataset.tab);
    }));

    function tile(label, value, sub, id) {
        return `<div class="tile"><div class="tile-label">${esc(label)}</div>` +
            `<div class="tile-value"${id ? ` id="${id}"` : ''}>${value}</div>` +
            `<div class="tile-sub"${id ? ` id="${id}-sub"` : ''}>${sub || ''}</div></div>`;
    }

    function renderModesStatic() {
        const s = state.summary;
        const net = s.network;
        const springs = net.springs === 'structure' ? ', structure-based springs' : '';
        let tiles = tile('Elastic network', `${net.nodes}`,
            `nodes · ${net.dof} degrees of freedom · ${METHOD_NAME[net.method]}, cutoff ${net.cutoff} Å${springs}`);
        tiles += tile('Selected modes', '', '', 't-sel');
        tiles += tile('Mean collectivity', '', '', 't-coll');
        if (s.structure.predicted) {
            tiles += tile('B-factor agreement', 'pLDDT', 'predicted model: no experimental B-factors');
        } else {
            const r = s.gnm && s.gnm.bfactor_corr;
            tiles += tile('B-factor agreement', r === null || r === undefined ? '–' : `r = ${fmt(r)}`,
                'GNM mean-square fluctuations vs experimental B-factors');
        }
        if (s.target) tiles += tile('Second conformation', '', '', 't-target');
        $('modes-tiles').innerHTML = tiles;
        $('cumulative-sub').textContent = s.modes.exact
            ? `Σ(1/λ) of the first M modes over all ${s.modes.n_total} non-zero modes. The marker follows the number of modes.`
            : `Large network: shares are relative to the first ${s.modes.eigvals.length} modes only.`;

        // ligands
        const ligs = s.ligands || [];
        $('ligands-block').style.display = ligs.length ? 'block' : 'none';
        const handling = { follow: 'moves with its pocket', network: 'network nodes', remove: 'removed' };
        $('ligands-table').innerHTML = '<thead><tr><th>Ligand / ion</th><th>Type</th><th class="num">Heavy atoms</th>' +
            '<th>Treatment</th><th>Pocket residues (≤ 4.5 Å)</th></tr></thead><tbody>' +
            ligs.map((l) => `<tr><td>${esc(l.label)}</td><td>${esc(l.kind)}</td><td class="num">${l.heavy_atoms}</td>` +
                `<td>${esc(handling[l.handling] || l.handling)}</td><td>${l.pocket.map((p) => `<span class="pill">${esc(p)}</span>`).join('') || '–'}</td></tr>`).join('') +
            '</tbody>';

        const notes = s.notes || [];
        $('notes-title').style.display = notes.length ? 'block' : 'none';
        $('notes').innerHTML = notes.map((n) => `<li>${esc(n)}</li>`).join('');

        // table view of the modes
        const m = s.modes;
        let cum = 0;
        let rows = '';
        for (let k = 0; k < m.eigvals.length; k++) {
            cum += m.variance_fraction[k];
            rows += `<tr data-mode="${k + 1}"><td class="num">${k + 1}</td><td class="num">${m.eigvals[k].toPrecision(4)}</td>` +
                `<td class="num">${pct(m.variance_fraction[k], 2)}</td><td class="num">${pct(cum)}</td>` +
                `<td class="num">${fmt(m.collectivity[k])}</td>` +
                (s.target ? `<td class="num">${fmt(s.target.overlap[k])}</td><td class="num">${fmt(s.target.cumulative_overlap[k])}</td>` : '') +
                `<td>${esc(m.top_residues[k])}</td></tr>`;
        }
        $('modes-table').innerHTML = '<thead><tr><th class="num">Mode</th><th class="num">Eigenvalue</th>' +
            '<th class="num">Share</th><th class="num">Cumulative</th><th class="num">Collectivity</th>' +
            (s.target ? '<th class="num">Overlap</th><th class="num">Cum. overlap</th>' : '') +
            '<th>Most mobile residues</th></tr></thead><tbody>' + rows + '</tbody>';
    }

    function updateKDependent() {
        const s = state.summary;
        if (!s) return;
        const st = modeStats(state.K);
        const K = st.K;
        const set = (id, v) => { const el = $(id); if (el) el.innerHTML = v; };
        set('t-sel', `1–${K}`);
        set('t-sel-sub', `${pct(st.share)} of the network's fluctuation${st.exact ? '' : ' (computed modes)'}`);
        set('t-coll', fmt(st.coll));
        set('t-coll-sub', st.localized.length ? `${st.localized.length} localised mode(s) below 0.15` : 'no localised modes among them');
        if (s.target) {
            set('t-target', fmt(st.cumOverlap));
            set('t-target-sub', `cumulative overlap · RMSD ${fmt(st.targetRmsd)} → ${fmt(st.residual)} Å with ${K} modes`);
        }
        $('modes-guide').innerHTML = guideHtml(st);
        document.querySelectorAll('#modes-table tbody tr').forEach((tr) =>
            tr.classList.toggle('selected', Number(tr.dataset.mode) <= K));
        updateModeCharts();
    }

    function guideHtml(st) {
        const s = state.summary;
        const items = [];
        const need = st.need;
        items.push(`The slowest mode alone carries ${pct(s.modes.variance_fraction[0])} of the network's fluctuation; ` +
            `${need[0.5] ? `${need[0.5]} modes reach 50%` : '50% is not reached within the computed modes'}` +
            `${need[0.75] ? ` and ${need[0.75]} reach 75%` : ''}. ` +
            'Because amplitudes scale with 1/√λ, the first modes dominate whatever number you choose.');
        if (st.localized.length) {
            items.push(`Mode${st.localized.length > 1 ? 's' : ''} ${st.localized.slice(0, 10).join(', ')} ` +
                `${st.localized.length > 1 ? 'are' : 'is'} localised (collectivity < 0.15). Check in the mode explorer which ` +
                'residues move: if it is a disordered tail or an artefactual loop, remove it (chain selection, pLDDT ' +
                'trimming) or try structure-based springs or RTB.');
        } else {
            items.push(`All ${st.K} selected modes are collective (collectivity ≥ 0.15).`);
        }
        if (s.target) {
            items.push(`Compared with ${esc(s.target.name)}: ${st.need70 ? `${st.need70} mode(s) reach a cumulative overlap of 0.7` : 'a cumulative overlap of 0.7 is not reached'}` +
                `${st.need90 ? `, ${st.need90} reach 0.9` : ''}. This is the most direct evidence for the number of modes to use.`);
        } else {
            items.push('If you have a second structure of this protein (e.g. ligand-bound), add it under ' +
                '<em>Compare with a second conformation</em> to see how many modes the real change needs.');
        }
        return `<strong>Choosing the number of modes (currently ${st.K})</strong><ul>${items.map((i) => `<li>${i}</li>`).join('')}</ul>`;
    }

    // ------------------------------------------------------------------ charts
    const css = (name) => getComputedStyle(document.body).getPropertyValue(name).trim();
    const PLOT_CONFIG = {
        displaylogo: false, responsive: true, displayModeBar: 'hover',
        modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines', 'zoomIn2d',
            'zoomOut2d', 'pan2d', 'zoom2d', 'hoverClosestCartesian', 'hoverCompareCartesian'],
        toImageButtonOptions: { format: 'png', scale: 2 },
    };

    function axis(title, extra = {}) {
        return Object.assign({
            title: { text: title, standoff: 8 },
            gridcolor: css('--chart-grid'),
            linecolor: css('--chart-axis'),
            zerolinecolor: css('--chart-axis'),
            tickcolor: css('--chart-axis'),
            ticks: 'outside',
            ticklen: 4,
            automargin: true,
        }, extra);
    }

    function layout(xTitle, yTitle, extra = {}) {
        const l = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, system-ui, sans-serif', size: 12, color: css('--chart-text') },
            margin: { l: 56, r: 14, t: 12, b: 44 },
            xaxis: axis(xTitle, extra.xaxis || {}),
            yaxis: axis(yTitle, extra.yaxis || {}),
            hoverlabel: { font: { family: 'Inter, system-ui, sans-serif' } },
            showlegend: false,
            bargap: 0.15,
        };
        const rest = Object.assign({}, extra);
        delete rest.xaxis;
        delete rest.yaxis;
        return Object.assign(l, rest);
    }

    function legendTop() {
        return { showlegend: true, legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom', bgcolor: 'rgba(0,0,0,0)' }, margin: { l: 56, r: 14, t: 34, b: 44 } };
    }

    function vline(x, label) {
        return {
            shapes: [{ type: 'line', xref: 'x', yref: 'paper', x0: x, x1: x, y0: 0, y1: 1, line: { color: css('--chart-ink'), width: 1 } }],
            annotations: label ? [{
                xref: 'x', yref: 'paper', x, y: 1, text: label, showarrow: false, xanchor: x > 50 ? 'right' : 'left',
                xshift: x > 50 ? -4 : 4, yanchor: 'top', font: { size: 11, color: css('--chart-ink') },
                bgcolor: document.body.classList.contains('light-mode') ? 'rgba(255,255,255,0.85)' : 'rgba(20,22,27,0.85)',
            }] : [],
        };
    }

    const plotted = new Set();
    function plot(id, traces, lay) {
        Plotly.react(id, traces, lay, PLOT_CONFIG);
        plotted.add(id);
    }

    function resizeCharts() {
        for (const id of plotted) {
            const el = $(id);
            if (el && el.offsetParent !== null) Plotly.Plots.resize(el);
        }
    }

    function selectedColors(n, K) {
        const c1 = css('--chart-1');
        const muted = css('--chart-muted');
        return Array.from({ length: n }, (_, k) => (k < K ? c1 : muted));
    }

    function residueAxis(labels) {
        return labels.map((l, i) => i + 1);
    }

    function renderCharts() {
        const s = state.summary;
        if (!s) return;
        const m = s.modes;
        const n = m.eigvals.length;
        const modesX = Array.from({ length: n }, (_, k) => k + 1);
        const cum = [];
        m.variance_fraction.reduce((a, v, k) => (cum[k] = a + v), 0);
        const hoverMode = modesX.map((k) => `Mode ${k}<br>eigenvalue ${m.eigvals[k - 1].toPrecision(4)}` +
            `<br>share ${pct(m.variance_fraction[k - 1], 2)}<br>collectivity ${fmt(m.collectivity[k - 1])}` +
            `<br>${esc(m.top_residues[k - 1])}`);

        plot('chart-cumulative', [{
            x: modesX, y: cum.map((v) => 100 * v), type: 'scatter', mode: 'lines+markers',
            line: { color: css('--chart-1'), width: 2 }, marker: { size: 5, color: css('--chart-1') },
            hovertemplate: 'first %{x} modes: %{y:.1f}%<extra></extra>', name: 'cumulative share',
        }], layout('Number of modes M', '% of fluctuation', { yaxis: { range: [0, 101] } }));

        plot('chart-share', [{
            x: modesX, y: m.variance_fraction.map((v) => 100 * v), type: 'bar',
            marker: { color: selectedColors(n, state.K), line: { width: 0 } },
            hovertext: hoverMode, hoverinfo: 'text', name: 'share',
        }], layout('Mode', '% of fluctuation'));

        plot('chart-collectivity', [{
            x: modesX, y: m.collectivity, type: 'bar',
            marker: { color: selectedColors(n, state.K) }, hovertext: hoverMode, hoverinfo: 'text', name: 'collectivity',
        }], layout('Mode', 'Collectivity', {
            yaxis: { range: [0, 1] },
            shapes: [{ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: 0.15, y1: 0.15, line: { color: css('--chart-text'), width: 1 } }],
            annotations: [{ xref: 'paper', yref: 'y', x: 1, y: 0.15, text: 'localised below', showarrow: false, xanchor: 'right', yanchor: 'bottom', font: { size: 11, color: css('--chart-text') } }],
        }));

        for (const id of ['chart-share', 'chart-collectivity']) {
            const el = $(id);
            if (!el._rocClick) {
                el.on('plotly_click', (ev) => { if (ev.points.length) setK(ev.points[0].x); });
                el._rocClick = true;
            }
        }

        renderMobility();

        const hasTarget = !!s.target;
        $('overlap-card').style.display = hasTarget ? 'block' : 'none';
        $('residual-card').style.display = hasTarget ? 'block' : 'none';
        if (hasTarget) {
            const t = s.target;
            $('overlap-sub').textContent = `${t.name}: ${t.matched_residues} matched C-alpha atoms, RMSD ${fmt(t.rmsd)} Å after superposition. ` +
                'Overlap = |cos| between a mode and the conformational change; 1 means the mode alone describes it.';
            plot('chart-overlap', [
                { x: modesX, y: t.overlap, type: 'bar', name: 'overlap of each mode', marker: { color: css('--chart-1') }, hovertemplate: 'mode %{x}: %{y:.2f}<extra></extra>' },
                { x: modesX, y: t.cumulative_overlap, type: 'scatter', mode: 'lines', name: 'cumulative overlap, modes 1–M', line: { color: css('--chart-2'), width: 2 }, hovertemplate: 'modes 1–%{x}: %{y:.2f}<extra></extra>' },
            ], layout('Mode / number of modes M', 'Overlap', Object.assign({ yaxis: { range: [0, 1.02] } }, legendTop())));
            plot('chart-residual', [{
                x: modesX, y: t.residual_rmsd, type: 'scatter', mode: 'lines+markers', line: { color: css('--chart-1'), width: 2 },
                marker: { size: 5, color: css('--chart-1') }, hovertemplate: '%{x} modes: %{y:.2f} Å<extra></extra>', name: 'reachable RMSD',
            }], layout('Number of modes M', 'RMSD to second conformation (Å)', { yaxis: { rangemode: 'tozero' } }));
        }
        updateModeCharts();
        renderModeProfile(Number($('mode-select').value) || 1);
        if (s.ensemble) renderEnsembleCharts();
    }

    function updateModeCharts() {
        const s = state.summary;
        if (!s) return;
        const n = s.modes.eigvals.length;
        const K = Math.min(state.K, n);
        const cumK = s.modes.variance_fraction.slice(0, K).reduce((a, b) => a + b, 0);
        const v = vline(K, `M = ${K} · ${pct(cumK)}`);
        Plotly.relayout('chart-cumulative', { shapes: v.shapes, annotations: v.annotations });
        Plotly.restyle('chart-share', { 'marker.color': [selectedColors(n, K)] });
        Plotly.restyle('chart-collectivity', { 'marker.color': [selectedColors(n, K)] });
        if (s.target) {
            const t = s.target;
            const vo = vline(K, `M = ${K} · ${fmt(t.cumulative_overlap[K - 1])}`);
            Plotly.relayout('chart-overlap', { shapes: vo.shapes, annotations: vo.annotations });
            const vr = vline(K, `${fmt(t.residual_rmsd[K - 1])} Å`);
            vr.shapes.push({ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: t.rmsd, y1: t.rmsd, line: { color: css('--chart-text'), width: 1 } });
            vr.annotations.push({ xref: 'paper', yref: 'y', x: 1, y: t.rmsd, text: 'input vs second conformation', showarrow: false, xanchor: 'right', yanchor: 'bottom', font: { size: 11, color: css('--chart-text') } });
            Plotly.relayout('chart-residual', { shapes: vr.shapes, annotations: vr.annotations });
        }
        renderMobility();
    }

    function anmMobility(K) {
        const m = state.summary.modes;
        const nres = m.profiles[0].length;
        const out = new Array(nres).fill(0);
        for (let k = 0; k < Math.min(K, m.profiles.length); k++) {
            const p = m.profiles[k];
            const w = 1 / m.eigvals[k];
            for (let i = 0; i < nres; i++) out[i] += p[i] * w;
        }
        return out;
    }

    function normalise(a) {
        const mean = a.reduce((x, y) => x + y, 0) / a.length;
        return a.map((v) => (mean > 0 ? v / mean : v));
    }

    function renderMobility() {
        const s = state.summary;
        const res = s.residues;
        const x = residueAxis(res.labels);
        const traces = [];
        const predicted = s.structure.predicted;
        if (!predicted && res.bfactor) {
            traces.push({ x, y: normalise(res.bfactor), name: 'experimental B-factor', type: 'scatter', mode: 'lines', line: { color: css('--chart-3'), width: 2 }, text: res.labels, hovertemplate: '%{text}: %{y:.2f}<extra>B-factor</extra>' });
        }
        if (res.gnm_msf) {
            traces.push({ x, y: normalise(res.gnm_msf), name: 'GNM, all modes', type: 'scatter', mode: 'lines', line: { color: css('--chart-2'), width: 2 }, text: res.labels, hovertemplate: '%{text}: %{y:.2f}<extra>GNM</extra>' });
        }
        traces.push({ x, y: normalise(anmMobility(state.K)), name: `selected modes 1–${state.K}`, type: 'scatter', mode: 'lines', line: { color: css('--chart-1'), width: 2 }, text: res.labels, hovertemplate: '%{text}: %{y:.2f}<extra>modes 1–' + state.K + '</extra>' });
        $('mobility-title').textContent = predicted ? 'Residue mobility (predicted model)' : 'Residue mobility vs B-factors';
        $('mobility-sub').textContent = (predicted
            ? 'pLDDT is a confidence score, not a B-factor, so no experimental comparison is shown. '
            : 'Each profile divided by its mean. ') +
            'Peaks in the selected-modes profile are what the ensemble will move most.';
        plot('chart-mobility', traces, layout('Residue (index)', 'Relative mobility', legendTop()));
    }

    // ------------------------------------------------------------------ mode explorer
    function fillModeSelect() {
        const m = state.summary.modes;
        const sel = $('mode-select');
        sel.innerHTML = m.eigvals.slice(0, Math.min(50, m.eigvals.length)).map((_, k) =>
            `<option value="${k + 1}">Mode ${k + 1} · ${pct(m.variance_fraction[k])} · κ ${fmt(m.collectivity[k])}</option>`).join('');
    }

    $('mode-select').addEventListener('change', () => showMode(Number($('mode-select').value)));
    $('mode-play').addEventListener('click', () => {
        if (!state.modeViewer) return;
        state.modePlaying = !state.modePlaying;
        if (state.modePlaying) state.modeViewer.animate({ loop: 'forward', interval: 90 });
        else state.modeViewer.stopAnimate();
        $('mode-play').textContent = state.modePlaying ? 'Pause' : 'Play';
    });

    function renderModeProfile(k) {
        const s = state.summary;
        const p = s.modes.profiles[k - 1];
        if (!p) return;
        const amp = p.map(Math.sqrt);
        const max = Math.max(...amp) || 1;
        $('mode-profile-title').textContent = `Displacement profile of mode ${k}`;
        $('mode-profile-sub').textContent = `Relative C-alpha amplitude (1 = most mobile residue). Collectivity ${fmt(s.modes.collectivity[k - 1])}; ` +
            `most mobile: ${s.modes.top_residues[k - 1]}.`;
        plot('chart-mode-profile', [{
            x: residueAxis(s.residues.labels), y: amp.map((v) => v / max), type: 'scatter', mode: 'lines',
            line: { color: css('--chart-1'), width: 2 }, fill: 'tozeroy', fillcolor: 'rgba(57,135,229,0.12)',
            text: s.residues.labels, hovertemplate: '%{text}: %{y:.2f}<extra></extra>', name: `mode ${k}`,
        }], layout('Residue (index)', 'Relative amplitude', { yaxis: { range: [0, 1.05] } }));
    }

    // sequential blue ramp (dataviz reference palette), low -> high
    const RAMP = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];

    function rampColor(t) {
        const light = document.body.classList.contains('light-mode');
        const ramp = light ? RAMP : RAMP.slice().reverse();
        t = Math.max(0, Math.min(1, t)) * (ramp.length - 1);
        const i = Math.min(ramp.length - 2, Math.floor(t));
        const f = t - i;
        const a = parseInt(ramp[i].slice(1), 16), b = parseInt(ramp[i + 1].slice(1), 16);
        const ch = (sh) => Math.round(((a >> sh) & 255) * (1 - f) + ((b >> sh) & 255) * f);
        return (ch(16) << 16) | (ch(8) << 8) | ch(0);
    }

    function residueValueMap(chains, nums, values) {
        const map = new Map();
        const max = Math.max(...values) || 1;
        values.forEach((v, i) => map.set(`${chains[i]}:${nums[i]}`, v / max));
        return map;
    }

    function styleViewer(v, valueMap) {
        v._valueMap = valueMap;
        recolorViewer(v);
    }

    function recolorViewer(v) {
        const map = v._valueMap;
        const mode = v === state.ensViewer ? $('ens-color').value : 'value';
        let cartoon;
        if (mode === 'chain') cartoon = { colorscheme: 'chain' };
        else if (mode === 'spectrum') cartoon = { color: 'spectrum' };
        else cartoon = { colorfunc: (atom) => rampColor(map ? (map.get(`${atom.chain}:${atom.resi}`) ?? 0) : 0.5) };
        v.setStyle({}, { cartoon: Object.assign({ arrows: true }, cartoon) });
        v.setStyle({ hetflag: true }, { stick: { radius: 0.22, colorscheme: 'orangeCarbon' }, sphere: { scale: 0.25, colorscheme: 'orangeCarbon' } });
        v.setStyle({ resn: ['HOH', 'WAT'] }, {});
    }

    async function showMode(k) {
        if (!state.summary) return;
        $('mode-select').value = String(k);
        renderModeProfile(k);
        const msg = $('mode-viewer-msg');
        msg.textContent = `Loading mode ${k}…`;
        msg.style.display = 'flex';
        try {
            const r = await fetch(`/api/jobs/${state.resultJob}/mode/${k}`);
            if (!r.ok) throw new Error((await r.json()).error || 'Mode not available');
            const txt = await r.text();
            if (!state.modeViewer) {
                state.modeViewer = $3Dmol.createViewer($('mode-viewer'), { backgroundColor: viewerBackground(), antialias: true });
            }
            const v = state.modeViewer;
            v.stopAnimate();
            v.clear();
            v.addModelsAsFrames(txt, 'pdb');
            const s = state.summary;
            const amp = s.modes.profiles[k - 1].map(Math.sqrt);
            styleViewer(v, residueValueMap(s.residues.chains, s.residues.resnums, amp));
            v.zoomTo();
            v.render();
            if (state.modePlaying) v.animate({ loop: 'forward', interval: 90 });
            msg.style.display = 'none';
        } catch (e) {
            msg.textContent = `Could not load mode ${k}: ${e.message}`;
        }
    }

    // ------------------------------------------------------------------ ensemble
    function renderEnsembleStatic() {
        const s = state.summary;
        const e = s.ensemble;
        const q = e.qc;
        const rm = e.rmsd;
        const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;
        const median = (a) => { const b = a.slice().sort((x, y) => x - y); return b.length ? b[Math.floor(b.length / 2)] : 0; };
        const medPoly = median(q.clashes_polymer);
        const medLig = median(q.clashes_ligand);
        let tiles = tile('Models', `${e.n_models}`, `${METHOD_NAME[e.method]} · ${e.method === 'clustenm' ? 'ClustENM' : (e.sampling === 'traverse' ? 'mode traversal' : 'random combination')} · seed ${e.seed}`);
        tiles += tile('C-alpha RMSD to input', `${fmt(mean(rm))} Å`, `range ${fmt(Math.min(...rm))}–${fmt(Math.max(...rm))} Å`);
        tiles += tile('Radius of gyration', `${fmt(mean(e.rg))} Å`, `input ${fmt(e.rg_ref)} Å · range ${fmt(Math.min(...e.rg))}–${fmt(Math.max(...e.rg))}`);
        if (q.peptide_bad_fraction !== null && q.peptide_bad_fraction !== undefined) {
            tiles += tile('Peptide bonds off by > 0.5 Å', pct(q.peptide_bad_fraction, 2), `mean |Δ| ${fmt(q.peptide_mean, 3)} Å · worst ${fmt(q.peptide_max)} Å`);
        }
        tiles += tile('New steric overlaps', `${medPoly}`, `median per model, heavy atoms < ${q.clash_distance} Å${q.clashes_ligand.some((v) => v > 0) ? ` · ligands: ${medLig}` : ''}`);
        tiles += tile('Representatives', `${e.medoids.length}`, `models ${e.medoids.join(', ')}`);
        $('ens-tiles').innerHTML = tiles;

        // quality guidance
        const items = [];
        const bad = q.peptide_bad_fraction || 0;
        if (e.method === 'clustenm') {
            items.push('ClustENM conformers were energy-minimised with a force field: bond geometry and contacts are physically refined.');
        } else {
            items.push('Every residue (and ligand) was moved as a rigid body, so bond lengths and angles inside residues are exactly those of the input; only the links between residues can stretch.');
            if (bad < 0.005) items.push(`Peptide bonds: ${pct(bad, 2)} deviate by more than 0.5 Å, good local geometry.`);
            else items.push(`Peptide bonds: ${pct(bad, 2)} deviate by more than 0.5 Å${q.worst_bonds && q.worst_bonds.length ? ` (worst: ${q.worst_bonds.slice(0, 4).map((w) => `${esc(w[0])} ${fmt(w[1])} Å`).join(', ')})` : ''}. Consider a smaller RMSD, structure-based springs or RTB.`);
            if (medPoly > 5) items.push(`A median of ${medPoly} new heavy-atom overlaps per model: side chains collide at this amplitude. Lower the RMSD, use RTB (respects packing), or minimise the models (e.g. ClustENM, OpenMM, Rosetta relax) before docking or MD.`);
            else items.push(`Steric overlaps: median ${medPoly} per model, the packing is preserved.`);
        }
        if (e.sampling === 'random' && e.method !== 'clustenm') items.push('The ensemble RMSF should follow the network prediction for the selected modes; large differences point to residues damped by the rigid-residue reconstruction.');
        $('qc-guide').innerHTML = `<strong>Quality of the ensemble</strong><ul>${items.map((i) => `<li>${i}</li>`).join('')}</ul>`;

        // representatives
        const sizes = e.cluster_sizes;
        $('reps-table').innerHTML = '<thead><tr><th class="num">Cluster</th><th class="num">Models</th><th class="num">Representative</th><th class="num">RMSD to input (Å)</th></tr></thead><tbody>' +
            e.medoids.map((mdl, i) => `<tr><td class="num">${i + 1}</td><td class="num">${sizes[i]}</td><td class="num">model ${mdl}</td><td class="num">${fmt(rm[mdl - 1])}</td></tr>`).join('') + '</tbody>';

        const ligs = (s.ligands || []).filter((l) => l.shift_mean !== undefined);
        $('ens-ligands-block').style.display = ligs.length ? 'block' : 'none';
        $('ens-ligands-table').innerHTML = '<thead><tr><th>Ligand / ion</th><th>Treatment</th><th class="num">Mean displacement (Å)</th><th class="num">Max displacement (Å)</th></tr></thead><tbody>' +
            ligs.map((l) => `<tr><td>${esc(l.label)}</td><td>${esc(l.handling === 'network' ? 'network nodes' : 'moves with its pocket')}</td><td class="num">${fmt(l.shift_mean)}</td><td class="num">${fmt(l.shift_max)}</td></tr>`).join('') + '</tbody>';
        $('clash-sub').textContent = `Heavy-atom pairs closer than ${q.clash_distance} Å that were not in contact in the input.`;
    }

    function renderEnsembleCharts() {
        const s = state.summary;
        const e = s.ensemble;
        const models = e.rmsd.map((_, i) => i + 1);
        const hover = models.map((i) => `model ${i}` + (e.labels ? ` (mode ${e.labels[i - 1][0]}, ${fmt(e.labels[i - 1][1])} Å)` : '') + `<br>cluster ${e.clusters[i - 1]}`);
        const rmsdShapes = [];
        const rmsdAnn = [];
        if (e.sampling === 'random' && e.method !== 'clustenm') {
            rmsdShapes.push({ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: s.parameters.rmsd, y1: s.parameters.rmsd, line: { color: css('--chart-text'), width: 1 } });
            rmsdAnn.push({ xref: 'paper', yref: 'y', x: 1, y: s.parameters.rmsd, text: 'requested average', showarrow: false, xanchor: 'right', yanchor: 'bottom', font: { size: 11, color: css('--chart-text') } });
        }
        plot('chart-rmsd', [{
            x: models, y: e.rmsd, type: 'scatter', mode: models.length > 150 ? 'markers' : 'lines+markers',
            line: { color: css('--chart-1'), width: 1 }, marker: { size: 8, color: css('--chart-1'), line: { width: 2, color: 'rgba(0,0,0,0)' } },
            hovertext: hover, hovertemplate: '%{hovertext}<br>%{y:.2f} Å<extra></extra>', name: 'RMSD',
        }], layout('Model', 'C-alpha RMSD (Å)', { yaxis: { rangemode: 'tozero' }, shapes: rmsdShapes, annotations: rmsdAnn }));

        plot('chart-rg', [{
            x: models, y: e.rg, type: 'scatter', mode: 'markers', marker: { size: 8, color: css('--chart-1') },
            hovertext: hover, hovertemplate: '%{hovertext}<br>%{y:.2f} Å<extra></extra>', name: 'Rg',
        }], layout('Model', 'Radius of gyration (Å)', {
            shapes: [{ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: e.rg_ref, y1: e.rg_ref, line: { color: css('--chart-text'), width: 1 } }],
            annotations: [{ xref: 'paper', yref: 'y', x: 1, y: e.rg_ref, text: 'input', showarrow: false, xanchor: 'right', yanchor: 'bottom', font: { size: 11, color: css('--chart-text') } }],
        }));

        const labels = e.residue_labels;
        const xr = residueAxis(labels);
        const rmsfTraces = [{ x: xr, y: e.rmsf, type: 'scatter', mode: 'lines', name: 'ensemble RMSF', line: { color: css('--chart-1'), width: 2 }, text: labels, hovertemplate: '%{text}: %{y:.2f} Å<extra>ensemble</extra>' }];
        if (e.pred_rmsf && e.pred_rmsf.length === labels.length) {
            rmsfTraces.push({ x: xr, y: e.pred_rmsf, type: 'scatter', mode: 'lines', name: 'expected from the network', line: { color: css('--chart-2'), width: 2 }, text: labels, hovertemplate: '%{text}: %{y:.2f} Å<extra>network</extra>' });
        }
        plot('chart-rmsf', rmsfTraces, layout('Residue (index)', 'RMSF (Å)', Object.assign({ yaxis: { rangemode: 'tozero' } }, rmsfTraces.length > 1 ? legendTop() : {})));

        const q = e.qc;
        const clashTraces = [{ x: models, y: q.clashes_polymer, type: 'bar', name: 'macromolecule', marker: { color: css('--chart-1') }, hovertemplate: 'model %{x}: %{y}<extra>macromolecule</extra>' }];
        const ligAny = q.clashes_ligand.some((v) => v > 0);
        if (ligAny) clashTraces.push({ x: models, y: q.clashes_ligand, type: 'bar', name: 'involving ligands', marker: { color: css('--chart-2') }, hovertemplate: 'model %{x}: %{y}<extra>ligands</extra>' });
        plot('chart-clashes', clashTraces, layout('Model', 'Overlapping atom pairs', Object.assign({ barmode: 'stack', yaxis: { rangemode: 'tozero' } }, ligAny ? legendTop() : {})));

        $('pairwise-card').style.display = e.pairwise ? 'block' : 'none';
        if (e.pairwise) {
            const light = document.body.classList.contains('light-mode');
            const ramp = light ? RAMP : RAMP.slice().reverse();
            plot('chart-pairwise', [{
                z: e.pairwise, x: models, y: models, type: 'heatmap',
                colorscale: ramp.map((c, i) => [i / (ramp.length - 1), c]),
                colorbar: { title: { text: 'Å', side: 'top' }, thickness: 10, outlinewidth: 0, tickfont: { color: css('--chart-text') } },
                hovertemplate: 'models %{x} and %{y}: %{z:.2f} Å<extra></extra>',
            }], layout('Model', 'Model', {
                xaxis: { constrain: 'domain' },
                yaxis: { autorange: 'reversed', scaleanchor: 'x', constrain: 'domain' },
                margin: { l: 56, r: 10, t: 12, b: 44 },
            }));
        }

        $('rama-card').style.display = e.rama ? 'block' : 'none';
        if (e.rama) {
            const ens = e.rama.ensemble;
            const ref = e.rama.reference;
            plot('chart-rama', [
                { x: ens.map((p) => p[0]), y: ens.map((p) => p[1]), type: 'scattergl', mode: 'markers', name: 'ensemble models', marker: { size: 4, color: css('--chart-1'), opacity: 0.45 }, hovertemplate: 'φ %{x:.0f}°, ψ %{y:.0f}°<extra>ensemble</extra>' },
                { x: ref.map((p) => p[0]), y: ref.map((p) => p[1]), type: 'scattergl', mode: 'markers', name: 'input structure', marker: { size: 7, symbol: 'circle-open', color: css('--chart-2'), line: { width: 1.5, color: css('--chart-2') } }, hovertemplate: 'φ %{x:.0f}°, ψ %{y:.0f}°<extra>input</extra>' },
            ], layout('φ (°)', 'ψ (°)', Object.assign({
                xaxis: { range: [-180, 180], dtick: 90, constrain: 'domain' },
                yaxis: { range: [-180, 180], dtick: 90, scaleanchor: 'x', constrain: 'domain' },
            }, legendTop())));
        }
    }

    async function loadEnsembleViewer() {
        const s = state.summary;
        const e = s.ensemble;
        stopEnsemble();
        const msg = $('ens-viewer-msg');
        msg.style.display = 'flex';
        msg.textContent = 'Loading the ensemble…';
        const useAll = e.n_models <= 200;
        const file = useAll ? s.files.ensemble : s.files.representatives;
        $('ens-viewer-note').textContent = useAll ? '' : `The ensemble has ${e.n_models} models; the viewer shows the ${e.medoids.length} representatives.`;
        try {
            const r = await fetch(fileUrl(file));
            const txt = await r.text();
            if (!state.ensViewer) {
                state.ensViewer = $3Dmol.createViewer($('ens-viewer'), { backgroundColor: viewerBackground(), antialias: true });
            }
            const v = state.ensViewer;
            v.clear();
            v.addModelsAsFrames(txt, 'pdb');
            const keys = e.residue_keys || [];
            styleViewer(v, residueValueMap(keys.map((k) => k[0]), keys.map((k) => k[1]), e.rmsf));
            v.zoomTo();
            state.ensFrames = useAll ? e.n_models : e.medoids.length;
            state.ensFrame = 0;
            $('ens-frame').max = state.ensFrames - 1;
            $('ens-frame').value = 0;
            await setEnsFrame(0);
            msg.style.display = 'none';
        } catch (err) {
            msg.textContent = `Could not load the ensemble: ${err.message}`;
        }
    }

    async function setEnsFrame(i) {
        const v = state.ensViewer;
        if (!v) return;
        state.ensFrame = i;
        await v.setFrame(i);
        v.render();
        $('ens-frame').value = i;
        const s = state.summary;
        const model = s.ensemble.n_models <= 200 ? i + 1 : s.ensemble.medoids[i];
        $('ens-frame-label').textContent = `model ${model} · ${fmt(s.ensemble.rmsd[model - 1])} Å`;
    }

    function stopEnsemble() {
        if (state.ensTimer) clearInterval(state.ensTimer);
        state.ensTimer = null;
        $('ens-play').textContent = 'Play';
    }

    $('ens-play').addEventListener('click', () => {
        if (!state.ensViewer || state.ensFrames < 2) return;
        if (state.ensTimer) { stopEnsemble(); return; }
        $('ens-play').textContent = 'Pause';
        state.ensTimer = setInterval(() => setEnsFrame((state.ensFrame + 1) % state.ensFrames), 180);
    });
    $('ens-frame').addEventListener('input', () => { stopEnsemble(); setEnsFrame(Number($('ens-frame').value)); });
    $('ens-color').addEventListener('change', () => { if (state.ensViewer) { recolorViewer(state.ensViewer); state.ensViewer.render(); } });

    // ------------------------------------------------------------------ files
    const FILE_DESC = {
        ensemble: 'All models, multi-model PDB (MODEL/ENDMDL).',
        representatives: 'Cluster medoids of the pairwise RMSD, one model per cluster.',
        reference: 'The structure exactly as modelled (waters removed, chains/ligand options applied).',
        clustenm_reference: 'ClustENM starting structure (hydrogens and missing atoms added by PDBFixer).',
        rmsf_pdb: 'Reference structure with the ensemble RMSF in the B-factor column (colour by B in PyMOL/ChimeraX).',
        nmd: 'Normal modes for VMD NMWiz (C-alpha).',
        modes_csv: 'Every mode: eigenvalue, share of fluctuation, collectivity, overlap with the second conformation.',
        residues_csv: 'Per residue: B-factor or pLDDT, GNM fluctuation, hinge flag, ensemble and expected RMSF.',
        models_csv: 'Per model: RMSD, radius of gyration, cluster, peptide-bond deviations, steric overlaps.',
        pairwise_csv: 'Pairwise C-alpha RMSD matrix between models.',
        zip: 'Everything above plus the run log and summary.json.',
    };

    function renderFiles() {
        const f = state.summary.files;
        const order = ['ensemble', 'representatives', 'reference', 'clustenm_reference', 'rmsf_pdb', 'nmd',
            'modes_csv', 'residues_csv', 'models_csv', 'pairwise_csv', 'zip'];
        $('files-list').innerHTML = order.filter((k) => f[k]).map((k) =>
            `<div class="file-row"><div><a href="${fileUrl(f[k], true)}">${esc(f[k])}</a>` +
            `<div class="file-desc">${esc(FILE_DESC[k])}</div></div></div>`).join('') +
            (state.summary.parameters.split_models && state.summary.ensemble
                ? '<div class="file-row"><div>One PDB file per model<div class="file-desc">Inside the .zip, in the models/ folder.</div></div></div>' : '');
    }

    window.addEventListener('resize', () => {
        for (const v of [state.modeViewer, state.ensViewer]) if (v) { v.resize(); v.render(); }
    });

    updateMethod();
    updateButtons();

    // Reopen the results of a job from the address bar (?job=<id>), e.g. after a reload.
    const params = new URLSearchParams(window.location.search);
    const linked = params.get('job');
    if (linked && /^[0-9a-f]{12}$/.test(linked)) {
        fetch(`/api/jobs/${linked}`).then((r) => (r.ok ? r.json() : null)).then((d) => {
            if (!d) return;
            if (d.state === 'running') {
                state.job = state.inputJob = linked;
                state.inputDirty = false;
                $('job-panel').style.display = 'block';
                state.pollTimer = setTimeout(poll, 200);
                updateButtons();
            } else if (d.state === 'done') {
                state.job = state.resultJob = state.inputJob = linked;
                state.inputDirty = false;
                updateButtons();
                showResults(d.summary);
            }
        });
    }
});
