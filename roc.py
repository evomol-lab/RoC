#!/usr/bin/env python3
"""The Room of Conformations (RoC): normal-mode based all-atom conformational ensembles.

Engine and command-line interface. The web interface (RoCGUI.py) runs this file in a
subprocess, one job per output directory:

    python roc.py structure.pdb --outdir results/                  # modes + ensemble
    python roc.py structure.pdb --outdir results/ --analyze-only   # normal modes only
    python roc.py --help

Workflow
  1. Read the structure with ProDy (PDB or mmCIF, first model, altloc A), remove waters,
     classify every residue as protein, nucleic acid, ligand or ion.
  2. Build an elastic network model (ANM on C-alpha atoms, or RTB on all heavy atoms with
     rigid residues / secondary-structure blocks) and compute its normal modes.
  3. Describe the modes (variance share, collectivity, GNM hinges and B-factor agreement,
     optional overlap with a second conformation) to support the choice of mode number.
  4. Sample conformations along the slowest modes (random Gaussian combination with
     amplitudes proportional to 1/sqrt(eigenvalue), or single-mode traversal), scaled to a
     chosen average C-alpha RMSD.
  5. Rebuild every residue (and ligand) as a rigid body fitted (Kabsch) onto its displaced
     network nodes, so bond lengths and angles inside residues and ligands are preserved.
  6. Analyse the ensemble with MDAnalysis (RMSD, RMSF, radius of gyration, pairwise RMSD,
     Ramachandran angles) and check its geometry (peptide bonds, new steric overlaps).
"""
import argparse
import json
import os
import sys
import time
import traceback
import warnings
import zipfile

warnings.filterwarnings('ignore', category=SyntaxWarning)   # ProDy docstrings on Python 3.12+
warnings.filterwarnings('ignore', category=DeprecationWarning)

import numpy as np
import scipy.linalg
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree
import prody as pr

__version__ = '2.0.0'

pr.confProDy(verbosity='none')

DENSE_LIMIT = 6000        # largest Hessian dimension diagonalised in full (all eigenvalues)
N_ANALYSIS_MODES = 100    # modes described in the analysis (and offered in the interface)
WINDOW = 2                # C-alpha neighbours on each side used to orient a residue (ANM)
PEPTIDE_BOND = 2.0        # C(i)-N(i+1) distance (A) that marks a covalent peptide link
PEPTIDE_TOL = 0.5         # peptide bonds stretched/compressed by more than this (A) are flagged
CLASH_DIST = 2.2          # heavy atoms closer than this (A) overlap...
BONDED_REF = 2.5          # ...unless they were already this close in the input (bonded)
BOND_GAMMA = 100.0        # RTB: stiffness factor of covalent links between rigid blocks
POCKET_CUTOFF = 5.0       # polymer heavy atoms within this distance (A) define a ligand pocket
GNM_CUTOFF = 7.3          # GNM contact cutoff (A) used for B-factors and hinges
DEFAULT_CUTOFF = {'anm': 15.0, 'rtb': 8.0, 'clustenm': 15.0}


class RocError(Exception):
    """An input or parameter problem, reported to the user without a traceback."""


def log(msg=''):
    print(msg, flush=True)


def progress(frac, msg):
    print(f'@@progress {frac:.2f} {msg}', flush=True)


def import_mdanalysis():
    """Import MDAnalysis even when an optional, broken h5py is installed.

    MDAnalysis imports h5py for the H5MD trajectory format (unused by RoC) and only
    tolerates an ImportError. An h5py built against another NumPy version raises
    ValueError instead ("numpy.dtype size changed"), which would stop the analysis; such an
    h5py is hidden so MDAnalysis treats H5MD as unavailable.
    """
    if 'MDAnalysis' not in sys.modules:
        try:
            import h5py  # noqa: F401
        except ImportError:
            pass
        except Exception as exc:
            for name in [m for m in sys.modules if m == 'h5py' or m.startswith('h5py.')]:
                del sys.modules[name]
            sys.modules['h5py'] = None
            log(f'Note: h5py is installed but cannot be imported ({type(exc).__name__}: {exc}). '
                'It is not needed by RoC and was ignored; to repair it, reinstall h5py for your '
                'NumPy version (e.g. pip install --upgrade h5py).')
    import MDAnalysis
    return MDAnalysis


# --------------------------------------------------------------------------- structure

def read_structure(path):
    """First model of a PDB/mmCIF file (altloc A), with author chain IDs."""
    lower = path.lower()
    if lower.endswith(('.cif', '.mmcif', '.cif.gz', '.mmcif.gz')):
        ag = pr.parseMMCIF(path, unite_chains=True, model=1)
    else:
        ag = pr.parsePDB(path, model=1)
    if ag is None or ag.numAtoms() == 0:
        raise RocError(f'No atoms could be read from {os.path.basename(path)}.')
    ag.setSegnames(np.array([''] * ag.numAtoms()))
    return ag


def looks_predicted(path, ag):
    """True for AlphaFold-like models, whose B-factor column holds pLDDT."""
    name = os.path.basename(path).upper()
    hint = name.startswith('AF-') or 'ALPHAFOLD' in name
    if not hint:
        try:
            if path.lower().endswith('.gz'):
                import gzip
                with gzip.open(path, 'rb') as fh:
                    head = fh.read(400000)
            else:
                with open(path, 'rb') as fh:
                    head = fh.read(400000)
            head = head.decode('latin-1').upper()
            hint = 'ALPHAFOLD' in head or 'PLDDT' in head
        except OSError:
            pass
    b = ag.getBetas()
    return bool(hint and b is not None and b.min() >= 0 and b.max() <= 100)


class Residue:
    __slots__ = ('kind', 'atoms', 'heavy', 'chain', 'resname', 'resnum', 'icode',
                 'hetero', 'names', 'ca', 'segment')

    @property
    def label(self):
        return f'{self.chain}:{self.resname}{self.resnum}{self.icode}'.strip(':')


class System:
    """An AtomGroup with every residue classified and polymer segments identified."""

    def __init__(self, ag):
        self.ag = ag
        n = ag.numAtoms()
        X = ag.getCoords()
        names = ag.getNames()
        water = ag.getFlags('water')
        ion = ag.getFlags('ion')
        protein = ag.getFlags('protein')
        nucleic = ag.getFlags('nucleic')
        hetero = ag.getFlags('hetatm') if ag.getFlags('hetatm') is not None else np.zeros(n, bool)
        hydrogen = ag.getFlags('hydrogen')
        self.heavy_mask = ~hydrogen
        self.residues = []
        for res in ag.getHierView().iterResidues():
            r = Residue()
            idx = res.getIndices()
            r.atoms = idx
            r.heavy = idx[~hydrogen[idx]]
            r.names = {nm: i for nm, i in zip(names[idx], idx)}
            r.chain = res.getChid()
            r.resname = res.getResname()
            r.resnum = int(res.getResnum())
            r.icode = res.getIcode() or ''
            r.hetero = bool(hetero[idx].all())
            r.ca = None
            r.segment = -1
            backbone = {'N', 'CA', 'C'} <= set(r.names)
            if water[idx].any():
                r.kind = 'water'
            elif (protein[idx].any() or backbone) and not ion[idx].any():
                r.kind = 'protein' if 'CA' in r.names else 'other'
                if r.kind == 'protein':
                    r.ca = r.names['CA']
            elif nucleic[idx].any():
                r.kind = 'nucleic'
            elif ion[idx].any() or len(r.heavy) == 1:
                r.kind = 'ion'
            else:
                r.kind = 'ligand'
            self.residues.append(r)

        # Polymer segments: consecutive residues of a chain joined by a covalent link.
        seg = -1
        prev = None
        for r in self.residues:
            if r.kind not in ('protein', 'nucleic'):
                continue
            linked = False
            if prev is not None and prev.kind == r.kind and prev.chain == r.chain:
                if r.kind == 'protein':
                    if 'C' in prev.names and 'N' in r.names:
                        linked = np.linalg.norm(X[prev.names['C']] - X[r.names['N']]) < PEPTIDE_BOND
                    else:
                        linked = np.linalg.norm(X[prev.ca] - X[r.ca]) < 4.3
                elif "O3'" in prev.names and 'P' in r.names:
                    linked = np.linalg.norm(X[prev.names["O3'"]] - X[r.names['P']]) < PEPTIDE_BOND
            if not linked:
                seg += 1
            r.segment = seg
            prev = r
        # A lone HETATM amino acid (free ligand such as GLU or a modified residue that is not
        # linked to anything) is a ligand, not a one-residue chain.
        counts = np.bincount([r.segment for r in self.residues if r.segment >= 0]) if seg >= 0 else []
        for r in self.residues:
            if r.segment >= 0 and counts[r.segment] == 1 and r.hetero:
                r.kind, r.ca, r.segment = 'ligand', None, -1

        self.protein = [r for r in self.residues if r.kind == 'protein']
        self.nucleic = [r for r in self.residues if r.kind == 'nucleic']
        self.polymer = [r for r in self.residues if r.kind in ('protein', 'nucleic')]
        self.ions = [r for r in self.residues if r.kind == 'ion']
        self.others = [r for r in self.residues if r.kind == 'other']
        ligand_res = [r for r in self.residues if r.kind == 'ligand']
        self.ligands = self._group_ligands(ligand_res, X)
        self.ca_idx = np.array([r.ca for r in self.protein], dtype=int)
        self.polymer_heavy = (np.concatenate([r.heavy for r in self.polymer])
                              if self.polymer else np.zeros(0, int))

    @staticmethod
    def _group_ligands(ligand_res, X):
        """Covalently linked ligand residues (glycans, multi-residue inhibitors) move together."""
        n = len(ligand_res)
        parent = list(range(n))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        if n > 1:
            atoms = np.concatenate([r.heavy for r in ligand_res])
            owner = np.concatenate([[k] * len(r.heavy) for k, r in enumerate(ligand_res)])
            for i, j in cKDTree(X[atoms]).query_pairs(PEPTIDE_BOND):
                a, b = find(owner[i]), find(owner[j])
                if a != b:
                    parent[a] = b
        groups = {}
        for k, r in enumerate(ligand_res):
            groups.setdefault(find(k), []).append(r)
        return list(groups.values())

    def counts(self):
        return {'protein_residues': len(self.protein), 'nucleotides': len(self.nucleic),
                'ligands': len(self.ligands), 'ions': len(self.ions),
                'chains': sorted({r.chain for r in self.polymer}),
                'atoms': int(self.ag.numAtoms())}


def ligand_label(group):
    first = group[0]
    names = '+'.join(r.resname for r in group)
    return f'{first.chain}:{names} {first.resnum}{first.icode}'.strip(':')


def prepare_system(path, params, notes):
    """Read, filter and classify the input structure."""
    ag = read_structure(path)
    predicted = looks_predicted(path, ag)
    n_water = int(ag.getFlags('water').sum())
    keep = ~ag.getFlags('water')
    if params.chains:
        wanted = [c.strip() for c in params.chains.split(',') if c.strip()]
        chain_mask = np.isin(ag.getChids(), wanted)
        if not chain_mask.any():
            raise RocError(f'None of the chains {", ".join(wanted)} exist in the structure '
                           f'(chains found: {", ".join(sorted(set(ag.getChids())))}).')
        keep &= chain_mask
    ag1 = ag[np.where(keep)[0]].copy()
    if n_water:
        notes.append(f'{n_water} water atoms were removed (crystallographic waters are not '
                     'part of the ensemble).')
    system = System(ag1)

    drop = np.zeros(ag1.numAtoms(), bool)
    if params.trim_plddt > 0:
        if not predicted:
            notes.append('pLDDT trimming was requested but the structure does not look like a '
                         'predicted model; the B-factor column was used anyway.')
        low = [r for r in system.protein if ag1.getBetas()[r.ca] < params.trim_plddt]
        for r in low:
            drop[r.atoms] = True
        if low:
            notes.append(f'{len(low)} residues with pLDDT < {params.trim_plddt:g} were removed '
                         'before building the network.')
    if params.ligands == 'remove' or params.method == 'clustenm':
        n_lig = len(system.ligands) + len(system.ions)
        for group in system.ligands:
            for r in group:
                drop[r.atoms] = True
        for r in system.ions:
            drop[r.atoms] = True
        if n_lig:
            why = ('ClustENM rebuilds the protein with a force field that has no parameters '
                   'for them' if params.method == 'clustenm' else 'as requested')
            notes.append(f'{n_lig} ligand(s)/ion(s) were removed ({why}).')
    if drop.any():
        system = System(ag1[np.where(~drop)[0]].copy())
    if len(system.protein) < 3:
        raise RocError('At least three protein residues with C-alpha atoms are needed.')
    if system.others:
        notes.append(f'{len(system.others)} incomplete residue(s) without a C-alpha atom will '
                     'follow their surroundings.')
    return system, predicted


# --------------------------------------------------------------------------- network models

def anm_hessian(X, cutoff, gamma=1.0, pair_gamma=None):
    """Sparse ANM Hessian, built with a k-d tree.

    With a uniform *gamma* it is identical to ProDy's ANM.buildHessian. *pair_gamma(i, j, r2)*
    may return one spring constant per node pair instead.
    """
    n = len(X)
    pairs = cKDTree(X).query_pairs(cutoff, output_type='ndarray')
    if len(pairs) == 0:
        raise RocError('The elastic network has no springs; increase the cutoff.')
    i, j = pairs[:, 0], pairs[:, 1]
    d = X[j] - X[i]
    r2 = (d * d).sum(1)
    g = np.full(len(i), float(gamma)) if pair_gamma is None else pair_gamma(i, j, r2)
    blk = -g[:, None, None] * d[:, :, None] * d[:, None, :] / r2[:, None, None]
    a, b = np.meshgrid(np.arange(3), np.arange(3), indexing='ij')
    rows = np.concatenate([(3 * i)[:, None, None] + a, (3 * j)[:, None, None] + a]).ravel()
    cols = np.concatenate([(3 * j)[:, None, None] + b, (3 * i)[:, None, None] + b]).ravel()
    vals = np.concatenate([blk, blk]).ravel()
    diag = np.zeros((n, 3, 3))
    np.add.at(diag, i, -blk)
    np.add.at(diag, j, -blk)
    dr = ((3 * np.arange(n))[:, None, None] + a).ravel()
    dc = ((3 * np.arange(n))[:, None, None] + b).ravel()
    rows = np.concatenate([rows, dr])
    cols = np.concatenate([cols, dc])
    vals = np.concatenate([vals, diag.ravel()])
    return sp.csr_matrix((vals, (rows, cols)), shape=(3 * n, 3 * n))


def rtb_projection(X, blocks):
    """Sparse orthonormal basis of the rigid-body motions of every block (Tama et al. 2000)."""
    rows, cols, vals = [], [], []
    col = 0
    order = np.argsort(blocks, kind='stable')
    bounds = np.flatnonzero(np.diff(blocks[order])) + 1
    for idx in np.split(order, bounds):
        x = X[idx] - X[idx].mean(0)
        m = np.zeros((3 * len(idx), 6))
        for k in range(3):
            m[k::3, k] = 1.0
        m[1::3, 3], m[2::3, 3] = -x[:, 2], x[:, 1]
        m[0::3, 4], m[2::3, 4] = x[:, 2], -x[:, 0]
        m[0::3, 5], m[1::3, 5] = -x[:, 1], x[:, 0]
        u, s, _ = np.linalg.svd(m, full_matrices=False)
        u = u[:, s > 1e-6 * s[0]]
        r = (3 * idx[:, None] + np.arange(3)).ravel()
        for j in range(u.shape[1]):
            rows.append(r)
            cols.append(np.full(len(r), col))
            vals.append(u[:, j])
            col += 1
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(3 * len(X), col))


def diagonalise(H, n_keep, title, use_prody_anm):
    """Lowest non-zero modes of H. Returns (ProDy model, all non-zero eigenvalues or None)."""
    dim = H.shape[0]
    if dim <= DENSE_LIMIT:
        if use_prody_anm:
            model = pr.ANM(title)
            model.setHessian(H.toarray())
            model.calcModes(n_modes=None)          # ProDy drops the zero modes itself
            evals = model.getEigvals()
            if len(evals) > n_keep:
                model = model[:n_keep]
                sub = pr.ANM(title)
                sub.setEigens(model.getArray(), model.getEigvals())
                model = sub
            return model, evals
        w, v = scipy.linalg.eigh(H.toarray() if sp.issparse(H) else H)
    else:
        k = min(n_keep + 12, dim - 2)
        w, v = eigsh(H.tocsc(), k=k, sigma=-1e-3, which='LM')
        order = np.argsort(w)
        w, v = w[order], v[:, order]
    nonzero = w > 1e-6 * max(1.0, abs(w).max())
    evals = w[nonzero]
    model = pr.NMA(title)
    model.setEigens(v[:, nonzero][:, :n_keep], evals[:n_keep])
    return model, (evals if dim <= DENSE_LIMIT else None)


def assign_dssp(system, ref_pdb, notes):
    """Secondary structure ('H', 'E' or '-') of every protein residue, from MDAnalysis DSSP."""
    try:
        mda = import_mdanalysis()
        from MDAnalysis.analysis.dssp import DSSP   # MDAnalysis >= 2.8
        u = mda.Universe(ref_pdb)
        prot = u.atoms[np.concatenate([r.atoms for r in system.protein])]
        codes = DSSP(prot).run().results.dssp[0]
        # DSSP keeps only the residues MDAnalysis recognises as protein (standard names), so
        # its codes are mapped back through their C-alpha atoms; the others count as coil.
        cas = prot.select_atoms('protein and name CA')
        if len(cas) != len(codes):
            raise ValueError('residue count mismatch')
        by_ca = dict(zip(cas.indices, codes))
        skipped = [r.label for r in system.protein if r.ca not in by_ca]
        if skipped:
            notes.append('DSSP does not handle non-standard residues; '
                         + ', '.join(skipped[:5]) + ('…' if len(skipped) > 5 else '')
                         + ' counted as coil.')
        return {id(r): str(by_ca.get(r.ca, '-')) for r in system.protein}
    except Exception as exc:  # incomplete backbone, unusual residues...
        notes.append(f'DSSP could not assign secondary structure ({exc}).')
        return None


def sse_blocks(system, ss):
    """Rigid blocks: each helix (>= 4 residues) or strand (>= 3), and every other residue."""
    blocks, bid, run = {}, 0, []

    def close_run():
        nonlocal bid
        code = ss.get(id(run[0]), '-') if run else '-'
        if run and code in 'HE' and len(run) >= (4 if code == 'H' else 3):
            for r in run:
                blocks[id(r)] = bid
            bid += 1
        else:
            for r in run:
                blocks[id(r)] = bid
                bid += 1

    current = None
    for r in system.polymer:
        code = ss.get(id(r), '-')
        key = (r.segment, code)
        if key != current or code not in 'HE':
            close_run()
            run, current = [], key
        run.append(r)
    close_run()
    return blocks, bid


def structure_based_gamma(system, ss, node_of_ca, n_nodes, gamma):
    """Per-pair spring constants from ProDy's GammaStructureBased (Lezon & Bahar 2010).

    Applies to C-alpha node pairs (connected residues x10, helix and sheet contacts x6);
    every other pair keeps *gamma*.
    """
    ca = system.ag[system.ca_idx].copy()
    ca.setSecstrs(np.array([ss.get(id(r), '-') for r in system.protein]))
    gsb = pr.GammaStructureBased(ca, gamma=gamma)
    ca_index = np.full(n_nodes, -1)
    ca_index[node_of_ca] = np.arange(len(node_of_ca))

    def pair_gamma(i, j, r2):
        g = np.full(len(i), float(gamma))
        ci, cj = ca_index[i], ca_index[j]
        both = np.where((ci >= 0) & (cj >= 0))[0]
        g[both] = [gsb.gamma(r2[k], ci[k], cj[k]) for k in both]
        return g
    return pair_gamma


class Network:
    """Elastic network nodes, rigid fitting groups and normal modes of a System.

    ANM: one node per C-alpha (plus P/C4'/C2 per nucleotide); each residue is later placed
    by a rigid fit to the displaced nodes of a window of +-2 residues.
    RTB: all heavy atoms are nodes and residues (or secondary-structure elements) are rigid
    blocks (Tama et al. 2000); springs between covalently bonded atoms of neighbouring
    blocks are BOND_GAMMA times stiffer, and each block is placed by a rigid fit to its own
    displaced atoms plus the backbone of the flanking residues.
    """

    def __init__(self, system, method, cutoff, gamma, ligand_mode, blocks_kind, springs,
                 ref_pdb, n_keep, notes):
        self.method = method
        X = system.ag.getCoords()
        network_ligands = ligand_mode == 'network'
        node_atoms, node_block = [], []
        groups1, groups2 = [], []   # (atoms to move, fit points): stage 1 nodes, stage 2 atoms
        self.block_info = None
        self.springs = 'uniform'

        def add_nodes(atoms, block):
            start = len(node_atoms)
            node_atoms.extend(atoms)
            node_block.extend([block] * len(atoms))
            return np.arange(start, start + len(atoms))

        segments = {}
        for r in system.polymer:
            segments.setdefault(r.segment, []).append(r)
        ss = None
        if (method == 'rtb' and blocks_kind == 'sse') or (method == 'anm' and springs == 'structure'):
            ss = assign_dssp(system, ref_pdb, notes)

        if method == 'rtb':
            if blocks_kind == 'sse' and ss is not None:
                bmap, n_blocks = sse_blocks(system, ss)
                self.block_info = f'{n_blocks} rigid blocks (secondary-structure elements and loop residues)'
            else:
                bmap = {id(r): k for k, r in enumerate(system.polymer)}
                n_blocks = len(system.polymer)
                self.block_info = f'{n_blocks} rigid blocks (one per residue)'
            members = {}
            for r in system.polymer:
                members.setdefault(bmap[id(r)], []).append(r)
            heavy_node = {}
            block_res = []
            for b, residues in members.items():
                heavy = np.concatenate([r.heavy for r in residues])
                for a, n in zip(heavy, add_nodes(list(heavy), b)):
                    heavy_node[a] = n
                block_res.append(residues)
            for residues in block_res:
                pts = [heavy_node[a] for r in residues for a in r.heavy]
                seg = segments[residues[0].segment]
                first, last = seg.index(residues[0]), seg.index(residues[-1])
                for nb in (first - 1, last + 1):
                    if 0 <= nb < len(seg):
                        pts += [heavy_node[seg[nb].names[n]] for n in ('N', 'CA', 'C', 'P', "C4'")
                                if n in seg[nb].names and seg[nb].names[n] in heavy_node]
                groups1.append((np.concatenate([r.atoms for r in residues]), np.array(pts)))
            next_block = n_blocks
        else:
            res_nodes = {}
            for r in system.polymer:
                if r.kind == 'protein':
                    sel = [r.ca]
                else:
                    sel = [r.names[a] for a in ('P', "C4'", 'C2') if a in r.names] or list(r.heavy[:1])
                res_nodes[id(r)] = add_nodes(sel, len(node_block))
            for seg in segments.values():
                width = WINDOW if seg[0].kind == 'protein' else 1
                for p, r in enumerate(seg):
                    window = seg[max(0, p - width):p + width + 1]
                    if sum(len(res_nodes[id(w)]) for w in window) < 3:
                        window = seg[max(0, p - width - 1):p + width + 2]
                    groups1.append((r.atoms, np.concatenate([res_nodes[id(w)] for w in window])))
            next_block = len(node_block)

        # node index of each protein C-alpha, in residue order
        pos = {a: n for n, a in enumerate(node_atoms)}
        self.ca_nodes = np.array([pos[r.ca] for r in system.protein], int)

        pocket_tree = cKDTree(X[system.polymer_heavy])
        for group in system.ligands + [[r] for r in system.ions]:
            atoms = np.concatenate([r.atoms for r in group])
            heavy = np.concatenate([r.heavy for r in group])
            if network_ligands:
                pts = add_nodes(list(heavy), next_block)
                next_block += 1
                if len(pts) < 3 and group[0].kind != 'ion':
                    # too few atoms to orient the ligand: add the nearest C-alpha nodes
                    dist = np.linalg.norm(X[system.ca_idx] - X[heavy].mean(0), axis=1)
                    pts = np.concatenate([pts, self.ca_nodes[np.argsort(dist)[:4]]])
                groups1.append((atoms, pts))
            else:
                groups2.append((atoms, pocket_atoms(pocket_tree, system, X[heavy])))
        for r in system.others:
            groups2.append((r.atoms, pocket_atoms(pocket_tree, system, X[r.heavy])))

        self.node_atoms = np.array(node_atoms, int)
        self.node_block = np.array(node_block, int)
        self.n_nodes = len(self.node_atoms)
        self.X_nodes = X[self.node_atoms]

        t0 = time.time()
        if method == 'rtb':
            blk = self.node_block

            def pair_gamma(i, j, r2):
                bonded = (r2 < BONDED_REF ** 2) & (blk[i] != blk[j])
                return np.where(bonded, gamma * BOND_GAMMA, gamma)
            H = anm_hessian(self.X_nodes, cutoff, gamma, pair_gamma)
            P = rtb_projection(self.X_nodes, blk)
            Hb = P.T @ (H @ P)
            self.dof = Hb.shape[0]
            model, evals = diagonalise(Hb.toarray() if self.dof <= DENSE_LIMIT else Hb.tocsr(),
                                       n_keep, 'RTB', use_prody_anm=False)
            self.model = pr.NMA('RTB')
            self.model.setEigens(P @ model.getArray(), model.getEigvals())
        else:
            pair_gamma = None
            if springs == 'structure' and ss is not None:
                pair_gamma = structure_based_gamma(system, ss, self.ca_nodes, self.n_nodes, gamma)
                self.springs = 'structure'
            H = anm_hessian(self.X_nodes, cutoff, gamma, pair_gamma)
            self.dof = H.shape[0]
            self.model, evals = diagonalise(H, n_keep, 'ANM', use_prody_anm=True)
        self.all_evals = evals
        self.seconds = time.time() - t0
        self.n_zero = (self.dof - len(evals)) if evals is not None else None
        if self.n_zero is not None and self.n_zero > 6:
            notes.append(f'The network has {self.n_zero} zero-frequency modes instead of 6: some '
                         'parts are not connected to the rest. Consider a larger cutoff.')

        self.stage1 = Rebuilder(self.X_nodes, X, groups1)
        self.stage2 = Rebuilder(X, X, groups2) if groups2 else None
        moved = np.zeros(len(X), bool)
        for move, _ in groups1 + groups2:
            moved[move] = True
        if not moved.all():
            raise RuntimeError(f'{(~moved).sum()} atoms are not assigned to a rigid group.')

    @property
    def V(self):
        return self.model.getArray()

    @property
    def evals(self):
        return self.model.getEigvals()

    def build(self, displacements, X_ref):
        """All-atom models for node displacements of shape (n_models, n_nodes, 3)."""
        out = np.empty((len(displacements),) + X_ref.shape)
        for i, d in enumerate(displacements):
            out[i] = rebuild(self.stage1, self.stage2, X_ref, self.X_nodes + d)
        return out


def pocket_atoms(tree, system, Xlig):
    """Polymer heavy atoms around a ligand (at least 6; the 12 nearest otherwise)."""
    near = set()
    for hits in tree.query_ball_point(Xlig, POCKET_CUTOFF):
        near.update(hits)
    if len(near) < 6:
        _, nearest = tree.query(Xlig.mean(0), k=min(12, tree.n))
        near = set(np.atleast_1d(nearest).tolist())
    return system.polymer_heavy[np.array(sorted(near), int)]


class Rebuilder:
    """Moves groups of atoms as rigid bodies that best fit (Kabsch) displaced fit points.

    Each residue (or block, or ligand) keeps its internal geometry: it is rotated and
    translated as a whole onto the new positions of its fit points.
    """

    def __init__(self, P_ref, X_ref, groups):
        self.G = len(groups)
        self.i_fit = np.concatenate([np.asarray(f, int) for _, f in groups])
        self.g_fit = np.concatenate([np.full(len(f), g) for g, (_, f) in enumerate(groups)])
        self.i_mov = np.concatenate([np.asarray(m, int) for m, _ in groups])
        self.g_mov = np.concatenate([np.full(len(m), g) for g, (m, _) in enumerate(groups)])
        self.starts = np.r_[0, np.cumsum([len(f) for _, f in groups])[:-1]]
        self.cnt = np.array([len(f) for _, f in groups], float)
        P = P_ref[self.i_fit]
        self.Pc = np.add.reduceat(P, self.starts, axis=0) / self.cnt[:, None]
        self.Pd = P - self.Pc[self.g_fit]
        cov = np.add.reduceat(self.Pd[:, :, None] * self.Pd[:, None, :], self.starts, axis=0)
        sv = np.linalg.svd(cov, compute_uv=False)
        # fewer than 3 points, or (nearly) collinear points: translation only
        self.rigid = (self.cnt >= 3) & (sv[:, 1] > 1e-3 * np.maximum(sv[:, 0], 1e-12))
        self.local = X_ref[self.i_mov] - self.Pc[self.g_mov]

    def __call__(self, Q, out):
        """Write the moved atoms into *out* given displaced fit points *Q*."""
        Qf = Q[self.i_fit]
        Qc = np.add.reduceat(Qf, self.starts, axis=0) / self.cnt[:, None]
        Hm = np.add.reduceat(self.Pd[:, :, None] * (Qf - Qc[self.g_fit])[:, None, :],
                             self.starts, axis=0)
        U, _, Vt = np.linalg.svd(Hm)
        V = np.transpose(Vt, (0, 2, 1))
        d = np.sign(np.linalg.det(np.matmul(V, np.transpose(U, (0, 2, 1)))))
        V[:, :, 2] *= d[:, None]
        R = np.matmul(V, np.transpose(U, (0, 2, 1)))
        R[~self.rigid] = np.eye(3)
        out[self.i_mov] = np.einsum('nij,nj->ni', R[self.g_mov], self.local) + Qc[self.g_mov]
        return out

    def arrays(self, prefix):
        return {f'{prefix}_{k}': getattr(self, k) for k in
                ('i_fit', 'g_fit', 'i_mov', 'g_mov', 'starts', 'cnt', 'Pc', 'Pd', 'rigid', 'local')}

    @classmethod
    def from_arrays(cls, data, prefix):
        obj = cls.__new__(cls)
        for k in ('i_fit', 'g_fit', 'i_mov', 'g_mov', 'starts', 'cnt', 'Pc', 'Pd', 'rigid', 'local'):
            setattr(obj, k, data[f'{prefix}_{k}'])
        obj.G = len(obj.cnt)
        return obj


def rebuild(net_stage1, net_stage2, X_ref, node_coords):
    """All-atom coordinates for one set of displaced network nodes."""
    out = X_ref.copy()
    net_stage1(node_coords, out)
    if net_stage2 is not None:
        net_stage2(out.copy(), out)
    return out


# --------------------------------------------------------------------------- mode analysis

def mode_profiles(V, ca_nodes):
    """Squared C-alpha displacement of every residue in every mode (n_res x n_modes)."""
    n_nodes = V.shape[0] // 3
    Vr = V.reshape(n_nodes, 3, -1)[ca_nodes]
    return (Vr ** 2).sum(1)


def collectivity(sq):
    """Brueschweiler (1995) collectivity of each mode from per-residue squared amplitudes."""
    p = sq / sq.sum(0, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        ent = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=0)
    return np.exp(ent) / sq.shape[0]


def gnm_analysis(system, predicted, notes):
    """GNM on C-alpha atoms: predicted mean-square fluctuations and hinge residues."""
    ca = system.ag.getCoords()[system.ca_idx]
    gnm = pr.GNM('GNM')
    try:
        gnm.buildKirchhoff(ca, cutoff=GNM_CUTOFF)
        gnm.calcModes(n_modes=None if len(ca) <= 3000 else 200)
    except Exception as exc:
        notes.append(f'GNM analysis failed: {exc}')
        return None
    msf = pr.calcSqFlucts(gnm)
    betas = system.ag.getBetas()[system.ca_idx]
    corr = None
    if not predicted and np.ptp(betas) > 0:
        corr = float(np.corrcoef(msf, betas)[0, 1])
    hinges = {}
    for k in range(min(2, gnm.numModes())):
        hinges[k + 1] = [int(h) for h in pr.calcHinges(gnm[k])]
    return {'msf': msf, 'bfactor_corr': corr, 'hinges': hinges}


def target_analysis(system, network, target_path, n_modes, notes):
    """Overlap of the modes with the deformation from the input to a second conformation."""
    try:
        tgt = read_structure(target_path)
    except Exception as exc:
        notes.append(f'The second conformation could not be read: {exc}')
        return None
    ref_ag = system.ag
    tgt_ca = tgt.select('protein and name CA')
    if tgt_ca is None:
        notes.append('The second conformation has no protein C-alpha atoms.')
        return None
    ca_pos = {a: k for k, a in enumerate(system.ca_idx)}
    try:
        matches = pr.matchChains(ref_ag.select('protein'), tgt.select('protein'),
                                 subset='calpha', seqid=60, overlap=60)
    except Exception as exc:
        notes.append(f'Chain matching with the second conformation failed: {exc}')
        return None
    if not matches:
        notes.append('No chain of the second conformation matches the input (>=60% identity).')
        return None
    matches = sorted(matches, key=lambda m: (-(m[2] * m[3]),
                                             m[0][0].getChid() != m[1][0].getChid()))
    used_r, used_t, pairs = set(), set(), []
    for m in matches:
        rc, tc = m[0][0].getChid(), m[1][0].getChid()
        if rc in used_r or tc in used_t:
            continue
        used_r.add(rc)
        used_t.add(tc)
        for a, b in zip(m[0].getIndices(), m[1].getCoords()):
            if a in ca_pos:
                pairs.append((ca_pos[a], b))
    if len(pairs) < 10:
        notes.append('Too few matching residues with the second conformation.')
        return None
    res_k = np.array([p[0] for p in pairs])
    A = ref_ag.getCoords()[system.ca_idx][res_k]
    B = np.array([p[1] for p in pairs])
    # superpose the target on the input (Kabsch)
    Ac, Bc = A.mean(0), B.mean(0)
    U, _, Vt = np.linalg.svd((B - Bc).T @ (A - Ac))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    B_fit = (B - Bc) @ R.T + Ac
    defo = (B_fit - A).ravel()
    rmsd = float(np.sqrt((defo ** 2).sum() / len(A)))
    nodes = network.ca_nodes[res_k]
    V = network.V.reshape(network.n_nodes, 3, -1)[nodes].reshape(3 * len(nodes), -1)[:, :n_modes]
    norms = np.linalg.norm(V, axis=0)
    overlap = np.abs(defo @ V) / (np.linalg.norm(defo) * norms)
    Q, _ = np.linalg.qr(V)
    captured = np.cumsum((Q.T @ defo) ** 2) / (defo @ defo)
    residual = np.sqrt(np.maximum(defo @ defo * (1 - captured), 0) / len(A))
    chains = sorted(used_t)
    return {'name': os.path.basename(target_path), 'matched_residues': len(A), 'rmsd': rmsd,
            'chains': chains, 'overlap': overlap, 'cumulative_overlap': np.sqrt(captured),
            'captured': captured, 'residual_rmsd': residual}


# --------------------------------------------------------------------------- sampling

def sample_random(V, evals, ca_nodes, n_confs, rmsd, rng):
    """NMSim / ProDy sampleModes scheme: R = R0 + s * sum_k r_k u_k / sqrt(lambda_k).

    r_k ~ N(0, 1); s makes the average C-alpha RMSD of the ensemble equal to *rmsd*.
    """
    r = rng.standard_normal((n_confs, V.shape[1])) / np.sqrt(evals)
    D = (r @ V.T).reshape(n_confs, -1, 3)
    ca_rmsd = np.sqrt((D[:, ca_nodes] ** 2).sum(2).mean(1))
    scale = rmsd / ca_rmsd.mean()
    return D * scale, scale


def sample_traverse(V, ca_nodes, frames, rmsd):
    """Each mode on its own, from -rmsd to +rmsd (C-alpha RMSD at both ends)."""
    out, labels = [], []
    amps = np.linspace(-1, 1, frames)
    for k in range(V.shape[1]):
        u = V[:, k].reshape(-1, 3)
        a = rmsd / np.sqrt((u[ca_nodes] ** 2).sum(1).mean())
        for t in amps:
            out.append(u * a * t)
            labels.append((k + 1, float(t * rmsd)))
    return np.array(out), labels


# --------------------------------------------------------------------------- ClustENM

def run_clustenm(system, params, notes):
    """ProDy ClustENM: ANM sampling + OpenMM energy minimisation (and optional MD)."""
    try:
        import openmm  # noqa: F401
        import pdbfixer  # noqa: F401
    except ImportError:
        raise RocError('ClustENM needs OpenMM and PDBFixer: pip install openmm pdbfixer')
    pr.confProDy(verbosity='info')
    ens = pr.ClustENM(system.ag.getTitle())
    ens.setAtoms(system.ag.select('not hydrogen').copy())
    ens.run(cutoff=params.cutoff, gamma=params.gamma, n_modes=params.modes,
            n_confs=params.confs, rmsd=params.rmsd, n_gens=params.clustenm_gens,
            maxclust=params.confs, sim=params.clustenm_md, platform='CPU', parallel=False)
    pr.confProDy(verbosity='none')
    atoms = ens.getAtoms().copy()
    coords = ens.getCoordsets()
    if coords is None or len(coords) == 0:
        raise RocError('ClustENM did not return any conformer.')
    ref = ens.getCoords()
    atoms.setCoords(ref)
    notes.append('ClustENM rebuilt the protein with PDBFixer (hydrogens and missing atoms '
                 'added) and minimised every conformer with OpenMM (Amber99SB-ILDN, implicit '
                 'solvent).')
    return atoms, ref, coords


# --------------------------------------------------------------------------- ensemble analysis

def geometry_qc(system, X_ref, coordsets):
    """Peptide-bond deviations and new heavy-atom overlaps of every model.

    Returns one dict per model and the peptide bonds with the largest deviations.
    """
    ag = system.ag
    pep, pep_labels = [], []
    for a, b in zip(system.protein[:-1], system.protein[1:]):
        if a.segment == b.segment and 'C' in a.names and 'N' in b.names:
            pep.append((a.names['C'], b.names['N']))
            pep_labels.append(f'{a.label}–{b.resname}{b.resnum}{b.icode}')
    pep = np.array(pep, int).reshape(-1, 2)
    d0 = np.linalg.norm(X_ref[pep[:, 0]] - X_ref[pep[:, 1]], axis=1) if len(pep) else None
    worst = np.zeros(len(pep))
    heavy = np.where(system.heavy_mask)[0]
    lig = np.zeros(ag.numAtoms(), bool)
    for group in system.ligands:
        for r in group:
            lig[r.atoms] = True
    for r in system.ions:
        lig[r.atoms] = True
    rows = []
    for X in coordsets:
        row = {}
        if d0 is not None and len(d0):
            dev = np.abs(np.linalg.norm(X[pep[:, 0]] - X[pep[:, 1]], axis=1) - d0)
            row['peptide_mean'] = float(dev.mean())
            row['peptide_max'] = float(dev.max())
            row['peptide_bad'] = int((dev > PEPTIDE_TOL).sum())
            worst = np.maximum(worst, dev)
        pairs = cKDTree(X[heavy]).query_pairs(CLASH_DIST, output_type='ndarray')
        n_poly = n_lig = 0
        if len(pairs):
            a, b = heavy[pairs[:, 0]], heavy[pairs[:, 1]]
            new = np.linalg.norm(X_ref[a] - X_ref[b], axis=1) >= BONDED_REF
            a, b = a[new], b[new]
            involves = lig[a] | lig[b]
            n_lig = int(involves.sum())
            n_poly = int((~involves).sum())
        row['clashes_polymer'] = n_poly
        row['clashes_ligand'] = n_lig
        rows.append(row)
    order = np.argsort(worst)[::-1][:8]
    worst_bonds = [(pep_labels[k], float(worst[k])) for k in order if worst[k] > PEPTIDE_TOL]
    return rows, worst_bonds, len(pep)


def mdanalysis_ensemble(ref_pdb, ens_pdb, ca_idx, protein_idx, n_clusters):
    """RMSD, RMSF, radius of gyration, pairwise RMSD and backbone dihedrals (MDAnalysis)."""
    mda = import_mdanalysis()
    from MDAnalysis.analysis import align, rms
    from MDAnalysis.analysis.diffusionmap import DistanceMatrix
    from MDAnalysis.coordinates.memory import MemoryReader
    from functools import partial
    warnings.filterwarnings('ignore')

    ref = mda.Universe(ref_pdb)
    u = mda.Universe(ref_pdb, ens_pdb)
    out = {}
    # One sequential pass over the multi-model PDB: random access to its frames (as in a
    # pairwise matrix) re-parses the file and takes minutes for a few hundred models.
    prot, ca = u.atoms[protein_idx], u.atoms[ca_idx]
    ca_xyz, rg = [], []
    for _ in u.trajectory:
        ca_xyz.append(ca.positions.copy())
        rg.append(prot.radius_of_gyration())
    out['rg'] = np.array(rg)
    out['rg_ref'] = float(ref.atoms[protein_idx].radius_of_gyration())

    # C-alpha atoms only, in memory: RMSD to the input, pairwise RMSD and RMSF
    ref_ca = mda.Merge(ref.atoms[ca_idx])
    u_ca = mda.Merge(ref.atoms[ca_idx])
    u_ca.load_new(np.array(ca_xyz), format=MemoryReader)
    out['rmsd'] = rms.RMSD(u_ca, ref_ca).run().results.rmsd[:, 2]
    n = u_ca.trajectory.n_frames
    if n >= 2:
        dm = DistanceMatrix(u_ca, metric=partial(rms.rmsd, center=True, superposition=True)).run()
        D = dm.results.dist_matrix
        out['pairwise'] = D
        out['clusters'], out['medoids'] = cluster_models(D, n_clusters)
    else:
        out['pairwise'] = np.zeros((1, 1))
        out['clusters'], out['medoids'] = np.zeros(1, int), [0]
    align.AlignTraj(u_ca, ref_ca, in_memory=True).run()
    out['rmsf'] = rms.RMSF(u_ca.atoms).run().results.rmsf

    try:
        from MDAnalysis.analysis.dihedrals import Ramachandran
        # Ramachandran accepts only residues in MDAnalysis' own 'protein' selection
        # (standard residue names): modified residues such as CCS or MSE are left out.
        ref_rama = Ramachandran(ref.atoms[protein_idx].select_atoms('protein')).run()
        rama = Ramachandran(u.atoms[protein_idx].select_atoms('protein')).run()
        out['rama_ref'] = ref_rama.results.angles[0]
        out['rama'] = rama.results.angles
        out['rama_error'] = None
    except Exception as exc:
        out['rama_ref'] = out['rama'] = None
        out['rama_error'] = f'Ramachandran angles could not be computed ({exc}).'
    return out


def cluster_models(D, k):
    """Average-linkage clustering of the pairwise RMSD matrix; medoid of every cluster."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    n = len(D)
    k = max(1, min(k, n))
    if k == 1 or n < 3:
        labels = np.ones(n, int)
    else:
        Z = linkage(squareform((D + D.T) / 2, checks=False), method='average')
        labels = fcluster(Z, t=k, criterion='maxclust')
    clusters = np.zeros(n, int)
    medoids = []
    order = sorted(set(labels), key=lambda c: -(labels == c).sum())
    for new, c in enumerate(order, 1):
        members = np.where(labels == c)[0]
        clusters[members] = new
        sub = D[np.ix_(members, members)]
        medoids.append(int(members[np.argmin(sub.sum(1))]))
    return clusters, medoids


# --------------------------------------------------------------------------- output helpers

def write_models(ag, coordsets, path):
    out = ag.copy()
    out.setCoords(np.asarray(coordsets))
    pr.writePDB(path, out)


def rnd(a, nd=4):
    if a is None:
        return None
    return np.round(np.asarray(a, float), nd).tolist()


def write_csv(path, header, rows):
    with open(path, 'w') as fh:
        fh.write(','.join(header) + '\n')
        for row in rows:
            fh.write(','.join('' if v is None else (f'{v:.5g}' if isinstance(v, float) else str(v))
                              for v in row) + '\n')


def mode_trajectory(jobdir, mode, frames=15, rmsd=2.0):
    """Multi-model PDB that animates one mode (used by the web interface)."""
    state = np.load(os.path.join(jobdir, 'state.npz'))
    V = state['V']
    if not 1 <= mode <= V.shape[1]:
        raise RocError(f'Mode {mode} is not available.')
    ref = pr.parsePDB(os.path.join(jobdir, str(state['reference'])))
    X_ref = ref.getCoords()
    s1 = Rebuilder.from_arrays(state, 's1')
    s2 = Rebuilder.from_arrays(state, 's2') if 's2_cnt' in state else None
    u = V[:, mode - 1].reshape(-1, 3)
    ca_nodes = state['ca_nodes']
    a = rmsd / np.sqrt((u[ca_nodes] ** 2).sum(1).mean())
    nodes = state['X_nodes']
    phases = np.sin(np.linspace(0, 2 * np.pi, frames, endpoint=False))
    coords = [rebuild(s1, s2, X_ref, nodes + u * a * t) for t in phases]
    path = os.path.join(jobdir, f'mode_{mode:03d}.pdb')
    write_models(ref, coords, path)
    return path


# --------------------------------------------------------------------------- main pipeline

def run(params):
    t_start = time.time()
    os.makedirs(params.outdir, exist_ok=True)
    notes = []
    name = params.name or os.path.splitext(os.path.basename(params.input))[0]
    name = ''.join(c if c.isalnum() or c in '-_.' else '_' for c in name).strip('.') or 'structure'
    if params.cutoff is None:
        params.cutoff = DEFAULT_CUTOFF[params.method]
    if params.seed is None:
        params.seed = int(np.random.SeedSequence().entropy % (2 ** 31))
    files = {}

    def path(suffix):
        return os.path.join(params.outdir, f'{name}_{suffix}')

    log(f'The Room of Conformations {__version__}  |  ProDy {pr.__version__}')
    log('─' * 72)
    progress(0.02, 'Reading structure')
    system, predicted = prepare_system(params.input, params, notes)
    c = system.counts()
    log(f'Structure: {os.path.basename(params.input)}  ({c["atoms"]} atoms after cleaning, '
        f'chains {", ".join(c["chains"])})')
    log(f'  protein residues: {c["protein_residues"]}   nucleotides: {c["nucleotides"]}   '
        f'ligands: {c["ligands"]}   ions: {c["ions"]}')
    if predicted:
        log('  predicted model detected: the B-factor column is read as pLDDT')

    ref_pdb = path('reference.pdb')
    pr.writePDB(ref_pdb, system.ag)
    files['reference'] = os.path.basename(ref_pdb)

    method_for_modes = 'anm' if params.method == 'clustenm' else params.method
    n_keep = max(params.modes, N_ANALYSIS_MODES)
    progress(0.08, 'Building the elastic network and computing normal modes')
    log(f'Elastic network: {method_for_modes.upper()}, cutoff {params.cutoff:g} Å, '
        f'gamma {params.gamma:g}, ligands: {params.ligands}')
    net = Network(system, method_for_modes, params.cutoff, params.gamma,
                  'remove' if params.method == 'clustenm' else params.ligands,
                  params.blocks, params.springs, ref_pdb, n_keep, notes)
    n_avail = net.model.numModes()
    log(f'  {net.n_nodes} nodes, {net.dof} degrees of freedom'
        + (f', {net.block_info}' if net.block_info else '')
        + (', structure-based springs' if net.springs == 'structure' else '')
        + f'; {n_avail} modes in {net.seconds:.1f} s')
    if params.modes > n_avail:
        notes.append(f'Only {n_avail} modes are available; using all of them.')
        params.modes = n_avail

    # ---- mode analysis
    progress(0.25, 'Analysing the normal modes')
    V_all, ev = net.V, net.evals
    sq = mode_profiles(V_all, net.ca_nodes)
    coll = collectivity(sq)
    inv = 1.0 / ev
    if net.all_evals is not None:
        total = (1.0 / net.all_evals).sum()
        exact = True
    else:
        total = inv.sum()
        exact = False
        notes.append(f'The network is large: variance shares are relative to the first '
                     f'{len(ev)} modes only.')
    frac = inv / total
    labels = [r.label for r in system.protein]
    top = []
    for k in range(len(ev)):
        order = np.argsort(sq[:, k])[::-1][:3]
        top.append(', '.join(labels[i] for i in order))
    gnm = gnm_analysis(system, predicted, notes)
    target = None
    if params.target:
        progress(0.30, 'Comparing with the second conformation')
        target = target_analysis(system, net, params.target, len(ev), notes)
        if target:
            log(f'Second conformation: {target["matched_residues"]} matched C-alpha atoms, '
                f'RMSD {target["rmsd"]:.2f} Å to the input')

    K = params.modes
    log(f'Selected modes: 1–{K}  →  {100 * frac[:K].sum():.1f}% of the network fluctuation'
        + ('' if exact else ' (of the computed modes)')
        + f'; mean collectivity {coll[:K].mean():.2f}')
    localized = [k + 1 for k in range(K) if coll[k] < 0.15]
    if localized:
        notes.append('Mode(s) ' + ', '.join(map(str, localized)) + ' are very localised '
                     '(collectivity < 0.15): they move only a few residues, often a flexible '
                     'terminus or loop. ' + ('Trimming low-pLDDT residues may help.'
                                             if predicted else ''))

    state = {'V': V_all.astype(np.float32), 'evals': ev, 'X_nodes': net.X_nodes,
             'ca_nodes': net.ca_nodes, 'reference': os.path.basename(ref_pdb)}
    state.update(net.stage1.arrays('s1'))
    if net.stage2 is not None:
        state.update(net.stage2.arrays('s2'))
    np.savez_compressed(os.path.join(params.outdir, 'state.npz'), **state)

    nmd = path('modes.nmd')
    try:
        nca = pr.NMA('modes')
        Vca = V_all.reshape(net.n_nodes, 3, -1)[net.ca_nodes].reshape(-1, V_all.shape[1])
        nca.setEigens(Vca[:, :max(K, 20)] / np.linalg.norm(Vca[:, :max(K, 20)], axis=0),
                      ev[:max(K, 20)])
        pr.writeNMD(nmd, nca, system.ag[system.ca_idx])
        files['nmd'] = os.path.basename(nmd)
    except Exception as exc:
        notes.append(f'NMD file not written: {exc}')

    mode_rows = []
    for k in range(len(ev)):
        mode_rows.append([k + 1, float(ev[k]), float(frac[k]), float(frac[:k + 1].sum()),
                          float(coll[k]),
                          float(target['overlap'][k]) if target else None,
                          float(target['cumulative_overlap'][k]) if target else None, top[k]])
    write_csv(path('modes.csv'), ['mode', 'eigenvalue', 'variance_fraction',
                                  'cumulative_variance', 'collectivity', 'overlap_target',
                                  'cumulative_overlap_target', 'top_residues'],
              [r[:-1] + [f'"{r[-1]}"'] for r in mode_rows])
    files['modes_csv'] = os.path.basename(path('modes.csv'))

    betas = system.ag.getBetas()[system.ca_idx]
    hinge_set = set()
    if gnm:
        for v in gnm['hinges'].values():
            hinge_set.update(v)
    summary = {
        'version': __version__,
        'action': 'analyze' if params.analyze_only else 'ensemble',
        'name': name,
        'parameters': {k: v for k, v in vars(params).items()
                       if k not in ('outdir',) and not callable(v)},
        'structure': dict(c, input=os.path.basename(params.input), predicted=predicted),
        'ligands': [],
        'network': {'method': method_for_modes, 'nodes': net.n_nodes, 'dof': net.dof,
                    'springs': net.springs,
                    'blocks': net.block_info, 'cutoff': params.cutoff, 'gamma': params.gamma,
                    'zero_modes': net.n_zero, 'seconds': round(net.seconds, 2)},
        'modes': {'n_total': int(len(net.all_evals)) if exact else None, 'exact': exact,
                  'eigvals': rnd(ev, 6), 'variance_fraction': rnd(frac, 6),
                  'collectivity': rnd(coll), 'top_residues': top,
                  'profiles': rnd(sq.T, 5)},
        'residues': {'labels': labels, 'chains': [r.chain for r in system.protein],
                     'resnums': [r.resnum for r in system.protein],
                     'bfactor': rnd(betas, 2),
                     'gnm_msf': rnd(gnm['msf'], 5) if gnm else None,
                     'hinge': [i in hinge_set for i in range(len(labels))]},
        'gnm': ({'bfactor_corr': gnm['bfactor_corr'],
                 'hinges': {k: [labels[i] for i in v] for k, v in gnm['hinges'].items()}}
                if gnm else None),
        'target': ({k: (rnd(v, 5) if isinstance(v, np.ndarray) else v)
                    for k, v in target.items()} if target else None),
        'notes': notes,
        'files': files,
    }
    X0 = system.ag.getCoords()
    pocket_tree = cKDTree(X0[system.polymer_heavy])
    owner = {}
    for r in system.polymer:
        for a in r.heavy:
            owner[a] = r
    for group in system.ligands + [[r] for r in system.ions]:
        heavy = np.concatenate([r.heavy for r in group])
        near = {}
        for hits in pocket_tree.query_ball_point(X0[heavy], 4.5):
            for i in hits:
                r = owner[system.polymer_heavy[i]]
                near[id(r)] = r
        pocket_labels = [r.label for r in system.polymer if id(r) in near]
        summary['ligands'].append({'label': ligand_label(group),
                                   'kind': 'ion' if group[0].kind == 'ion' else 'ligand',
                                   'heavy_atoms': int(len(heavy)),
                                   'handling': params.ligands,
                                   'pocket': pocket_labels})

    if params.analyze_only:
        finish(params, summary, system, betas, gnm, None, t_start, files, path)
        return summary

    # ---- ensemble
    rng = np.random.default_rng(params.seed)
    X_ref = system.ag.getCoords()
    model_labels = None
    pred_rmsf = None
    if params.method == 'clustenm':
        progress(0.35, 'Running ClustENM (ANM sampling + OpenMM minimisation)')
        log(f'ClustENM: {K} modes, {params.confs} conformers per generation, '
            f'{params.clustenm_gens} generation(s), MD: {"yes" if params.clustenm_md else "no"}')
        atoms_out, X_ref, coordsets = run_clustenm(system, params, notes)
        out_system = System(atoms_out)
        # ClustENM adds hydrogens and missing atoms: the ensemble gets its own reference
        ref_pdb = path('clustenm_reference.pdb')
        write_models(atoms_out, [X_ref], ref_pdb)
        files['clustenm_reference'] = os.path.basename(ref_pdb)
    else:
        progress(0.35, 'Sampling conformations along the normal modes')
        V, evK = V_all[:, :K], ev[:K]
        if params.sampling == 'traverse':
            D, model_labels = sample_traverse(V, net.ca_nodes, params.confs, params.rmsd)
            series = [np.where(np.array([m for m, _ in model_labels]) == k + 1)[0]
                      for k in range(K)]
            log(f'Traversing modes 1–{K}: {params.confs} frames each, ±{params.rmsd:g} Å '
                f'(C-alpha RMSD at the ends)')
        else:
            D, scale = sample_random(V, evK, net.ca_nodes, params.confs, params.rmsd, rng)
            series = [np.arange(len(D))]
            log(f'Random combination of modes 1–{K}: {params.confs} conformations, average '
                f'C-alpha RMSD {params.rmsd:g} Å (seed {params.seed})')
        progress(0.45, 'Rebuilding all-atom models')
        coordsets = net.build(D, X_ref)
        # Residues are placed by rigid fits to several nodes, which damps very localised
        # displacements slightly: rescale once so the rebuilt models reach the requested RMSD.
        ca = system.ca_idx
        achieved = np.sqrt(((coordsets[:, ca] - X_ref[ca]) ** 2).sum(2).mean(1))
        factors = np.ones(len(D))
        for idx in series:
            got = achieved[idx].mean() if params.sampling == 'random' else achieved[idx].max()
            if got > 0:
                factors[idx] = params.rmsd / got
        if np.abs(factors - 1).max() > 0.02:
            D = D * factors[:, None, None]
            coordsets = net.build(D, X_ref)
        if params.sampling == 'random':
            pred_rmsf = factors[0] * scale * np.sqrt((sq[:, :K] / evK).sum(1))
        out_system = system

    ens_pdb = path('ensemble.pdb')
    progress(0.6, 'Writing the ensemble')
    write_models(out_system.ag, coordsets, ens_pdb)
    files['ensemble'] = os.path.basename(ens_pdb)
    log(f'Ensemble: {len(coordsets)} models written to {os.path.basename(ens_pdb)}')

    progress(0.7, 'Checking geometry')
    qc, worst_bonds, n_pep = geometry_qc(out_system, X_ref, coordsets)
    for info, group in zip(summary['ligands'], out_system.ligands + [[r] for r in out_system.ions]):
        heavy = np.concatenate([r.heavy for r in group])
        shift = np.sqrt(((coordsets[:, heavy] - X_ref[heavy]) ** 2).sum(2).mean(1))
        info['shift_mean'] = round(float(shift.mean()), 3)
        info['shift_max'] = round(float(shift.max()), 3)
    progress(0.78, 'Analysing the ensemble with MDAnalysis')
    prot_idx = np.concatenate([r.atoms for r in out_system.protein])
    md = mdanalysis_ensemble(ref_pdb, ens_pdb, out_system.ca_idx, prot_idx, params.clusters)
    if md['rama_error']:
        notes.append(md['rama_error'])

    reps = path('representatives.pdb')
    write_models(out_system.ag, coordsets[md['medoids']], reps)
    files['representatives'] = os.path.basename(reps)
    if params.split_models:
        mdir = os.path.join(params.outdir, 'models')
        os.makedirs(mdir, exist_ok=True)
        for i, X in enumerate(coordsets):
            write_models(out_system.ag, [X], os.path.join(mdir, f'{name}_model_{i + 1:04d}.pdb'))
    rmsf_pdb = path('rmsf.pdb')
    out_ag = out_system.ag.copy()
    out_ag.setCoords(X_ref)
    per_atom = np.zeros(out_ag.numAtoms())
    for r, v in zip(out_system.protein, md['rmsf']):
        per_atom[r.atoms] = v
    out_ag.setBetas(per_atom)
    pr.writePDB(rmsf_pdb, out_ag)
    files['rmsf_pdb'] = os.path.basename(rmsf_pdb)

    n = len(coordsets)
    model_rows = []
    for i in range(n):
        q = qc[i]
        row = [i + 1, float(md['rmsd'][i]), float(md['rg'][i]), int(md['clusters'][i]),
               q.get('peptide_mean'), q.get('peptide_max'), q.get('peptide_bad'),
               q['clashes_polymer'], q['clashes_ligand']]
        if model_labels:
            row += [model_labels[i][0], model_labels[i][1]]
        model_rows.append(row)
    header = ['model', 'rmsd_ca', 'radius_of_gyration', 'cluster', 'peptide_bond_mean_dev',
              'peptide_bond_max_dev', f'peptide_bonds_off_by_{PEPTIDE_TOL}A',
              'new_overlaps_polymer', 'new_overlaps_ligand']
    if model_labels:
        header += ['mode', 'displacement_rmsd']
    write_csv(path('models.csv'), header, model_rows)
    files['models_csv'] = os.path.basename(path('models.csv'))
    D = md['pairwise']
    write_csv(path('pairwise_rmsd.csv'), ['model'] + [str(i + 1) for i in range(len(D))],
              [[i + 1] + [float(x) for x in D[i]] for i in range(len(D))])
    files['pairwise_csv'] = os.path.basename(path('pairwise_rmsd.csv'))

    pep_max = [q.get('peptide_max') for q in qc if q.get('peptide_max') is not None]
    pep_mean = [q.get('peptide_mean') for q in qc if q.get('peptide_mean') is not None]
    rama = None
    if md['rama'] is not None:
        r = md['rama']
        step = max(1, int(np.ceil(r.shape[0] * r.shape[1] / 20000)))
        flat = r.reshape(-1, 2)[::step]
        rama = {'reference': rnd(md['rama_ref'], 1), 'ensemble': rnd(flat, 1)}
    summary['ensemble'] = {
        'n_models': n, 'method': params.method, 'sampling': params.sampling,
        'seed': params.seed,
        'rmsd': rnd(md['rmsd'], 3), 'rg': rnd(md['rg'], 3), 'rg_ref': round(md['rg_ref'], 3),
        'rmsf': rnd(md['rmsf'], 3), 'pred_rmsf': rnd(pred_rmsf, 3),
        'residue_labels': [r.label for r in out_system.protein],
        'residue_keys': [[r.chain, r.resnum, r.icode] for r in out_system.protein],
        'clusters': md['clusters'].tolist(), 'medoids': [m + 1 for m in md['medoids']],
        'cluster_sizes': [int((md['clusters'] == c).sum()) for c in range(1, len(md['medoids']) + 1)],
        'pairwise': rnd(D, 2) if n <= 1000 else None,
        'qc': {'peptide_mean': float(np.mean(pep_mean)) if pep_mean else None,
               'peptide_max': float(np.max(pep_max)) if pep_max else None,
               'peptide_bad_fraction': (float(sum(q.get('peptide_bad', 0) for q in qc))
                                        / (n_pep * n) if n_pep else None),
               'peptide_tolerance': PEPTIDE_TOL, 'clash_distance': CLASH_DIST,
               'worst_bonds': worst_bonds,
               'clashes_polymer': [q['clashes_polymer'] for q in qc],
               'clashes_ligand': [q['clashes_ligand'] for q in qc]},
        'labels': [list(x) for x in model_labels] if model_labels else None,
        'rama': rama,
    }
    log(f'  C-alpha RMSD to input: mean {np.mean(md["rmsd"]):.2f} Å, '
        f'max {np.max(md["rmsd"]):.2f} Å')
    if pep_max:
        bad = summary['ensemble']['qc']['peptide_bad_fraction']
        log(f'  Peptide bonds: mean |Δ| {np.mean(pep_mean):.3f} Å, worst {np.max(pep_max):.2f} Å, '
            f'{100 * bad:.2f}% off by more than {PEPTIDE_TOL} Å')
    log(f'  New heavy-atom overlaps (< {CLASH_DIST} Å): median '
        f'{np.median([q["clashes_polymer"] for q in qc]):.0f} per model (polymer)'
        + (f', {np.median([q["clashes_ligand"] for q in qc]):.0f} (ligands)'
           if system.ligands or system.ions else ''))
    log(f'  {len(md["medoids"])} cluster representatives: models '
        + ', '.join(str(m + 1) for m in md['medoids']))
    same_residues = len(out_system.protein) == len(system.protein)
    finish(params, summary, system, betas, gnm, md['rmsf'] if same_residues else None,
           t_start, files, path, pred_rmsf)
    return summary


def finish(params, summary, system, betas, gnm, rmsf, t_start, files, path, pred_rmsf=None):
    rows = []
    for i, r in enumerate(system.protein):
        rows.append([r.chain, r.resname, r.resnum, r.icode, float(betas[i]),
                     float(gnm['msf'][i]) if gnm else None,
                     summary['residues']['hinge'][i],
                     float(rmsf[i]) if rmsf is not None and i < len(rmsf) else None,
                     float(pred_rmsf[i]) if pred_rmsf is not None else None])
    write_csv(path('residues.csv'), ['chain', 'resname', 'resnum', 'icode', 'bfactor',
                                     'gnm_msf', 'gnm_hinge', 'ensemble_rmsf',
                                     'enm_predicted_rmsf'], rows)
    files['residues_csv'] = os.path.basename(path('residues.csv'))
    summary['seconds'] = round(time.time() - t_start, 1)
    progress(0.95, 'Packing results')
    zpath = path('RoC_results.zip')
    with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED) as z:
        for f in sorted(set(files.values())):
            z.write(os.path.join(params.outdir, f), f)
        mdir = os.path.join(params.outdir, 'models')
        if os.path.isdir(mdir):
            for f in sorted(os.listdir(mdir)):
                z.write(os.path.join(mdir, f), f'models/{f}')
        z.writestr(f'{summary["name"]}_summary.json', json.dumps(summary, indent=1))
        log_file = os.path.join(params.outdir, 'run.log')
        if os.path.exists(log_file):
            z.write(log_file, f'{summary["name"]}_run.log')
    files['zip'] = os.path.basename(zpath)
    with open(os.path.join(params.outdir, 'summary.json'), 'w') as fh:
        json.dump(summary, fh)
    for note in summary['notes']:
        log(f'Note: {note}')
    log('─' * 72)
    log(f'Done in {summary["seconds"]:.1f} s.')
    progress(1.0, 'Done')


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog='roc.py', description='The Room of Conformations: all-atom conformational '
        'ensembles from elastic network normal modes (ProDy + MDAnalysis).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('input', help='structure file (.pdb, .cif, optionally .gz)')
    p.add_argument('--outdir', default='roc_results', help='output directory')
    p.add_argument('--name', help='prefix of the output files (default: input file name)')
    p.add_argument('--analyze-only', action='store_true',
                   help='compute and describe the normal modes, without an ensemble')
    p.add_argument('--target', help='second conformation of the same protein (e.g. ligand-bound '
                   'form): reports how well the modes describe the change')
    g = p.add_argument_group('model')
    g.add_argument('--method', choices=['anm', 'rtb', 'clustenm'], default='anm',
                   help='anm: C-alpha anisotropic network; rtb: rotation-translation of rigid '
                   'blocks on all heavy atoms; clustenm: ProDy ClustENM (needs OpenMM)')
    g.add_argument('--cutoff', type=float, help='spring cutoff in Å (default 15 for ANM, 8 for RTB)')
    g.add_argument('--gamma', type=float, default=1.0, help='spring constant')
    g.add_argument('--springs', choices=['uniform', 'structure'], default='uniform',
                   help='ANM spring constants: uniform (classic ANM) or structure-based '
                   '(ProDy GammaStructureBased: stiffer covalent, helix and sheet contacts)')
    g.add_argument('--blocks', choices=['residue', 'sse'], default='residue',
                   help='RTB rigid blocks: residues, or secondary-structure elements (DSSP)')
    g.add_argument('--ligands', choices=['follow', 'network', 'remove'], default='follow',
                   help='follow: ligands move rigidly with their binding pocket; network: '
                   'ligand heavy atoms are nodes of the elastic network; remove: drop them')
    g.add_argument('--chains', help='comma-separated chain IDs to keep (default: all)')
    g.add_argument('--trim-plddt', type=float, default=0,
                   help='remove residues with pLDDT below this value (predicted models)')
    g = p.add_argument_group('sampling')
    g.add_argument('--modes', type=int, default=15, help='number of slowest modes to combine')
    g.add_argument('--confs', type=int, default=50,
                   help='conformations (random) or frames per mode (traverse)')
    g.add_argument('--rmsd', type=float, default=1.0,
                   help='average C-alpha RMSD to the input (random) or RMSD at the ends of '
                   'each traversal (traverse), in Å')
    g.add_argument('--sampling', choices=['random', 'traverse'], default='random')
    g.add_argument('--seed', type=int, help='random seed (recorded in the results)')
    g.add_argument('--clusters', type=int, default=5, help='number of representative models')
    g.add_argument('--split-models', action='store_true', help='also write one PDB per model')
    g.add_argument('--clustenm-gens', type=int, default=1, help='ClustENM generations')
    g.add_argument('--clustenm-md', action='store_true',
                   help='ClustENM: short MD after minimisation (much slower)')
    args = p.parse_args(argv)
    if args.modes < 1:
        p.error('--modes must be at least 1')
    if not 1 <= args.confs <= 5000:
        p.error('--confs must be between 1 and 5000')
    if not 0 < args.rmsd <= 20:
        p.error('--rmsd must be between 0 and 20 Å')
    if args.sampling == 'traverse' and args.confs < 2:
        p.error('--confs must be at least 2 for traversal')
    return args


def main(argv=None):
    params = parse_args(argv)
    os.makedirs(params.outdir, exist_ok=True)
    status = os.path.join(params.outdir, 'status.json')
    try:
        run(params)
        state = {'state': 'done'}
        code = 0
    except RocError as exc:
        log(f'ERROR: {exc}')
        state = {'state': 'failed', 'error': str(exc)}
        code = 1
    except Exception as exc:
        log(traceback.format_exc())
        state = {'state': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
        code = 1
    with open(status, 'w') as fh:
        json.dump(state, fh)
    return code


if __name__ == '__main__':
    sys.exit(main())
