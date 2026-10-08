# The Room of Conformations: documentation

The Room of Conformations (RoC) generates all-atom conformational ensembles of proteins (and their
complexes) from the normal modes of elastic network models. It is a local, open replacement for
web servers such as NMSim (Krüger, Ahmed & Gohlke, *Nucleic Acids Res* 2012): it runs on your own
machine, keeps your structures private, and explains *why* the ensemble looks the way it does.

Contents

1. [Quick start](#quick-start)
2. [How RoC works](#how-roc-works)
3. [Elastic network models](#elastic-network-models)
4. [How many modes?](#how-many-modes)
5. [Amplitude and sampling schemes](#amplitude-and-sampling-schemes)
6. [All-atom reconstruction and quality checks](#all-atom-reconstruction-and-quality-checks)
7. [Ligands, ions, waters and other molecules](#ligands-ions-waters-and-other-molecules)
8. [Ensemble analysis](#ensemble-analysis)
9. [Output files](#output-files)
10. [Command line](#command-line)
11. [Comparison with NMSim and limitations](#comparison-with-nmsim-and-limitations)
12. [References](#references)

---

## Quick start

**Web interface**

```
pip install -r requirements.txt
python RoCGUI.py            # then open http://127.0.0.1:5050
```

1. **Structure**: upload a `.pdb` or `.cif` file (optionally gzipped) or fetch it from the RCSB PDB
   or AlphaFold DB. Optionally restrict the chains, trim low-pLDDT residues of predicted models, and
   add a *second conformation* of the same protein to compare with.
2. **Elastic network model**: choose ANM, RTB or ClustENM and how ligands are treated.
3. **Sampling**: choose the number of modes, the number of conformations and their average RMSD.
4. Click **Analyze modes** first (a few seconds): it computes the modes and shows the charts that
   help to choose the number of modes, without generating an ensemble. Adjust, then click
   **Generate ensemble**.

Results stay on the page; the address bar holds `?job=<id>`, so a reload reopens them. Every job is
stored in the `roc_jobs/` folder next to where the server was started (set `ROC_JOBS_DIR` to change
it, `ROC_JOB_RETENTION_DAYS` to delete old jobs automatically, `ROC_HOST`/`ROC_PORT` to change the
address).

**Command line**

```
python roc.py structure.pdb --outdir results --modes 10 --confs 100 --rmsd 1.0
python roc.py structure.pdb --outdir results --analyze-only --target holo.pdb
```

---

## How RoC works

1. **Read and clean the structure (ProDy).** First model only, alternate location A, author chain
   IDs for mmCIF. Waters are removed. Every residue is classified as protein, nucleic acid, ligand,
   ion or incomplete residue (see [Ligands](#ligands-ions-waters-and-other-molecules)).
2. **Build the elastic network and compute its normal modes (ProDy).** The Hessian is built with a
   k-d tree (identical to ProDy's `ANM.buildHessian`, much faster for large systems) and
   diagonalised in full up to 6000 degrees of freedom, so the share of every mode in the total
   fluctuation is exact. Larger systems use a sparse eigensolver for the slowest modes.
3. **Describe the modes** (variance share, collectivity, GNM B-factor agreement and hinges, overlap
   with a second conformation) so the number of modes is chosen on evidence, not by habit.
4. **Sample** conformations along the slowest modes.
5. **Rebuild all atoms**: every residue, ligand and ion is moved as a rigid body.
6. **Analyse the ensemble (MDAnalysis)** and check its geometry.

---

## Elastic network models

An elastic network model (ENM) replaces the protein by nodes (atoms or residues) connected by
harmonic springs of equal rest length to their distance in the input structure. Its normal modes
are the eigenvectors of the Hessian matrix; the eigenvalue λ of a mode is its stiffness. The six
zero-eigenvalue modes (rigid translation and rotation) are discarded, leaving 3*N* − 6 internal
modes for *N* nodes.

### ANM: anisotropic network model (default)

One node per C-alpha atom (plus P, C4' and C2 for each nucleotide), springs between nodes closer
than the **cutoff** (default 15 Å), spring constant **γ** (default 1) (Atilgan *et al.* 2001). γ
scales all eigenvalues equally and does not change the shape of the modes or of the ensemble.

**Spring constants**

- *Uniform*: every spring has the same γ: the classic ANM.
- *Structure-based*: ProDy's `GammaStructureBased` (Lezon & Bahar 2010): springs between
  covalently bonded residues are 10× stiffer, those between hydrogen-bonded residues of a helix or
  sheet 6× stiffer. Secondary structure comes from the MDAnalysis implementation of DSSP. In our
  tests this reduced the number of floppy, localised modes and the steric overlaps of the
  ensembles (e.g. trypsin 3PTB, HIV protease 1HSG) at no cost.

### RTB: rotations-translations of blocks

All heavy atoms are nodes (cutoff 8 Å by default) and groups of atoms move as **rigid blocks**
(Tama *et al.* 2000), either one block per residue or one block per helix (≥ 4 residues) and strand
(≥ 3 residues) with every loop residue on its own. The second choice is close in spirit to the
rigid-cluster normal mode analysis (RCNMA) used by NMSim. Springs between covalently bonded atoms
of neighbouring blocks (the peptide link) are 100× stiffer, which keeps the chain continuous.

Because RTB sees every atom, its modes respect side-chain packing. In our tests on compact
proteins (trypsin 3PTB, HIV protease 1HSG) at an average RMSD of 1.0 Å, C-alpha ANM ensembles had a
median of 10–16 new steric overlaps per model (7–10 with structure-based springs) and RTB ensembles
**none**. **RTB is the recommended method for ensembles meant for docking.** It costs a little more
(seconds for a few hundred residues) and its eigenvalues are not comparable with ANM eigenvalues.

### GNM: Gaussian network model (analysis only)

The isotropic GNM (Bahar *et al.* 1997; C-alpha atoms, cutoff 7.3 Å) gives no directions, so it
cannot generate conformations, but it predicts residue mean-square fluctuations very well. RoC uses
it to (i) report the correlation between predicted fluctuations and experimental B-factors, a quick
check that the network is a reasonable model of the structure, and (ii) list hinge residues of the
two slowest GNM modes (`calcHinges`).

### ClustENM (optional)

ProDy's ClustENM (Kurkcuoglu, Bahar & Doruker 2016) samples along ANM modes, clusters the
conformers and energy-minimises every representative with OpenMM (Amber99SB-ILDN with implicit
solvent; optional short MD), after PDBFixer has added hydrogens and missing atoms. It gives the
best stereochemistry, but takes minutes to hours and **removes ligands** (the force field has no
parameters for them). It needs `pip install openmm pdbfixer`; the option is disabled otherwise.

---

## How many modes?

### What a mode is

The modes are sorted by eigenvalue. The slowest modes (smallest λ) are the softest directions of
the network: large, collective motions such as hinge bending or the opening and closing of domains.
Faster modes are progressively more local. In thermal equilibrium the mean-square amplitude of a
mode is proportional to 1/λ, so a few slow modes carry a large part of the total fluctuation.

### What the number of modes changes, and what it does not

Every conformation is built as

  **R** = **R**₀ + *s* · Σₖ₌₁ᴹ *r*ₖ **u**ₖ / √λₖ,  with *r*ₖ ~ N(0, 1),

the scheme of ProDy's `sampleModes` and of NMSim's unbiased simulations: *M* is the number of
modes, **u**ₖ the mode vectors and *s* a single scale factor chosen so that the ensemble reaches the
requested average C-alpha RMSD. Therefore:

- **The amplitude is set by the RMSD, not by M.** Adding modes does not make the ensemble move
  more; it spreads the same deformation over more directions.
- **The first modes always dominate**, because of the 1/√λ weights. Going from 10 to 30 modes
  usually changes the ensemble much less than going from 1 to 10.
- **More modes = more local detail.** Higher modes add loop and side-chain-scale fluctuations and
  dilute the large collective motions.

### Rules of thumb

| Modes | What you get | Typical use |
|---|---|---|
| 1–3 | the dominant global motions only | studying one mechanism; with *Traverse*, a movie of each mode |
| 5–20 | the essential collective space | ensembles for docking or pocket analysis (default 15) |
| 30–100 | increasingly local "breathing" | approaching a thermal-like ensemble; diluted collective motion |

Ligand-induced conformational changes are typically well described by a handful of the slowest
modes (Tama & Sanejouand 2001; Bakan & Bahar 2009), which is why 5–20 modes is the usual range.

### Choosing on evidence: the charts

**Analyze modes** draws four diagnostics. The number of modes can be changed with the slider and
every chart, tile and the guidance box update immediately.

1. **Cumulative share of the fluctuation.** Σ(1/λ) of the first *M* modes divided by Σ(1/λ) of all
   modes: the fraction of the network's thermal motion described by your selection. A steep curve
   that flattens early means a few modes are enough.
2. **Collectivity** of each mode (Brüschweiler 1995): the effective fraction of residues that move,
   from 1/*N* (one residue) to 1 (all residues equally). Modes below 0.15 are flagged as
   *localised*: usually a disordered terminus, a long loop or an artefact (an isolated residue with
   few springs). Including them spends amplitude on wagging a tail. Inspect them in the **mode
   explorer**; remove the tail (chain selection or pLDDT trimming), or use structure-based springs
   or RTB.
3. **Residue mobility**: experimental B-factors (crystal structures), GNM prediction and the
   fluctuation profile of the selected modes, all divided by their mean. Peaks of the selected-modes
   profile are what the ensemble will move most.
4. **Overlap with a second conformation** (optional, and the most direct evidence). Upload another
   structure of the same protein (e.g. the ligand-bound form). RoC matches the chains (ProDy
   `matchChains`), superposes the C-alpha atoms and computes, for every mode, the overlap
   |cos θ| between the mode and the conformational change (Marques & Sanejouand 1995), the
   cumulative overlap of modes 1…*M* and the RMSD to the second conformation that remains after the
   best possible move along those modes. Choose *M* where the cumulative overlap levels off.

**Worked example: adenylate kinase** (open 4AKE chain A → closed 1AKE, C-alpha RMSD 7.1 Å, ANM with
default parameters). Mode 1 alone has an overlap of 0.80 with the closure; modes 1–5 reach a
cumulative overlap of 0.94 and could bring the structure to 2.4 Å of the closed form; modes 1–10
reach 0.97 (1.8 Å). Beyond 10 modes the gain is marginal: 10–15 modes are a sound choice here.

### Pitfalls

- **Predicted models (AlphaFold).** Low-confidence tails and linkers have few springs and dominate
  the slowest modes with localised motions. RoC detects AlphaFold models (pLDDT in the B-factor
  column) and suggests trimming, e.g. remove residues with pLDDT < 50 or 70.
- **Several chains.** In a complex the slowest modes are often rigid-body motions of the chains
  relative to each other. That is correct physics, but check that it is what you want; select a
  single chain otherwise. Do not compare asymmetric units whose chain packing differs (e.g. the
  dimeric asymmetric units of 1AKE and 4AKE): compare matching chains.
- **Large systems.** Above 6000 degrees of freedom only the slowest 100 modes are computed and the
  shares are relative to those modes.

---

## Amplitude and sampling schemes

**Random combination** (default) draws the coefficients of every conformation independently, as
above. The **average C-alpha RMSD** to the input sets the amplitude. 0.5–1.5 Å keeps local
geometry realistic for the C-alpha ANM; RTB tolerates larger values. Larger amplitudes explore
further, but the linear modes start to distort the structure: watch the quality report. Because
residues are rebuilt by rigid fits (next section), very localised displacements are slightly
damped; RoC rescales the displacements once so the rebuilt models reach the requested RMSD.

**Traverse each mode** moves along each of the first *M* modes on its own, from −RMSD to +RMSD in
the chosen number of frames: a movie of every mode, useful to understand what each mode does or to
build a mode-by-mode ensemble.

The **seed** of the random generator is recorded in the results; the same seed and parameters
give the same ensemble.

---

## All-atom reconstruction and quality checks

The normal modes move the network nodes. RoC then moves **every residue as a rigid body** (Kabsch
1976 superposition):

- **ANM**: the residue is placed by the best fit of its own C-alpha and its neighbours' C-alpha
  atoms (±2 residues in the same chain segment) onto their displaced positions. This captures the
  local rotation of the backbone, not only its translation.
- **RTB**: each block is placed by the best fit of its own heavy atoms plus the backbone of the
  flanking residues onto their displaced positions (this also removes the distortion of linearised
  rotations).
- **Hydrogens**, when present, move with their residue.

Bond lengths and angles **inside** residues and ligands are therefore exactly those of the input.
Only the links **between** residues can stretch. RoC reports for every model:

- **Peptide bonds**: the change of every C(i)–N(i+1) distance; the percentage of bonds off by more
  than 0.5 Å and the worst bonds are listed. At 1.0 Å average RMSD, fewer than 1% of the bonds
  exceeded 0.5 Å in all our tests (ANM and RTB), concentrated in a few exposed loops.
- **New steric overlaps**: heavy-atom pairs closer than 2.2 Å that were not in contact (≥ 2.5 Å)
  in the input, for the macromolecule and for ligands.

> Earlier versions of RoC (1.x) moved only the C-alpha atoms: every other atom stayed at its input
> position, so N–CA and CA–C bonds were stretched by up to several Å. Ensembles made with RoC 1.x
> should be regenerated.

The linear modes and rigid residues do not relax side chains. For docking or MD, minimising the
models is recommended when the report shows many overlaps (or use ClustENM).

---

## Ligands, ions, waters and other molecules

**Classification.** Every residue is classified with ProDy's atom flags and the structure itself:

| Class | Rule | Treatment |
|---|---|---|
| Protein | ProDy `protein` flag, or a residue with N, CA and C (modified residues such as MSE, SEP, TPO) | network nodes and rigid residues |
| Nucleic acid | ProDy `nucleic` flag | nodes P, C4', C2; rigid nucleotides |
| Water | ProDy `water` flag | always removed |
| Ion | ProDy `ion` flag or a single heavy atom | as ligands (below) |
| Ligand | anything else; covalently linked ligand residues (e.g. glycans) form one group | see options |
| Incomplete residue | protein-like but no C-alpha | moves with its surroundings |

A lone HETATM amino acid that is not bonded to the chain (a free amino-acid ligand) is a ligand.
Calcium ions named `CA` are recognised as ions, not as C-alpha atoms.

**Options**

- **Move with the binding site** (default). The ligand is not part of the network: the dynamics
  are those of the protein alone (apo-like). In every model the ligand is superposed rigidly onto
  its pocket, the polymer heavy atoms within 5 Å of it in the input. It keeps its internal
  geometry and follows its pocket, but nothing prevents the pocket from closing on it: check the
  ligand overlaps in the quality report.
- **Include in the elastic network.** Every ligand heavy atom and ion becomes a node. The ligand
  stiffens its binding site (holo-like dynamics) and moves with the modes; it is placed by a rigid
  fit to its own displaced nodes. Note that a ligand with many atoms adds many springs: this is the
  intended stiffening, but it weighs more than one residue.
- **Remove.** Only the macromolecule is modelled and written.

The results list every ligand and ion with its pocket residues (≤ 4.5 Å) and, after an ensemble,
its mean and maximum displacement.

---

## Ensemble analysis

The ensemble file is read back with MDAnalysis:

- **RMSD** of every model to the input (C-alpha atoms, after superposition).
- **RMSF** per residue after aligning all models on the input, compared with the RMSF expected from
  the network for the selected modes and amplitude.
- **Radius of gyration** of every model (protein atoms), against the input.
- **Pairwise RMSD** between models (`DistanceMatrix`), average-linkage clustering and the medoid
  of every cluster as **representative models**: a small, diverse subset for docking.
- **Ramachandran** φ/ψ angles of the ensemble and of the input.

---

## Output files

| File | Content |
|---|---|
| `*_ensemble.pdb` | all models (MODEL/ENDMDL) |
| `*_representatives.pdb` | cluster medoids |
| `models/*_model_NNNN.pdb` | one file per model (option, inside the zip) |
| `*_reference.pdb` | the structure exactly as modelled (waters removed, options applied) |
| `*_rmsf.pdb` | reference with the ensemble RMSF in the B-factor column |
| `*_modes.nmd` | normal modes for VMD's NMWiz plugin |
| `*_modes.csv` | eigenvalue, share, cumulative share, collectivity, overlap of every mode |
| `*_residues.csv` | per residue: B-factor/pLDDT, GNM fluctuation, hinge flag, ensemble and expected RMSF |
| `*_models.csv` | per model: RMSD, Rg, cluster, peptide-bond deviations, overlaps |
| `*_pairwise_rmsd.csv` | pairwise RMSD matrix |
| `*_summary.json`, `*_run.log` | every number shown in the interface, and the log |

---

## Command line

```
python roc.py INPUT [--outdir DIR] [--name PREFIX] [--analyze-only] [--target FILE]
              [--method {anm,rtb,clustenm}] [--cutoff Å] [--gamma G]
              [--springs {uniform,structure}] [--blocks {residue,sse}]
              [--ligands {follow,network,remove}] [--chains A,B] [--trim-plddt X]
              [--modes M] [--confs N] [--rmsd Å] [--sampling {random,traverse}]
              [--seed S] [--clusters K] [--split-models]
              [--clustenm-gens G] [--clustenm-md]
```

`python roc.py --help` describes every option and its default.

---

## Comparison with NMSim and limitations

| | NMSim web server | RoC |
|---|---|---|
| Network | rigid-cluster NMA (FIRST + RTB) | ANM, RTB (residue or secondary-structure blocks), GNM analysis |
| Sampling | iterative: modes recomputed after each step | single step from the input structure |
| Stereochemistry | constraint-based correction every step | rigid residues + quality report; ClustENM for force-field refinement |
| Guided simulations | targeted and radius-of-gyration guided | overlap analysis with a second conformation |
| Ligands | not modelled | follow the pocket, or part of the network |
| Where it runs | remote server, hours per job | your computer, seconds for most proteins |

Limitations to keep in mind:

- Normal modes are linear: very large amplitudes distort the structure. Large transitions (beyond
  2–3 Å) are better explored iteratively (NMSim, ClustENM with several generations) or with MD.
- Side chains are not repacked; models with overlaps should be minimised before docking or MD.
- The network sees only the coordinates you give it: crystal contacts, missing loops and membranes
  are ignored.

---

## References

- Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin O, Bahar I (2001) Anisotropy of
  fluctuation dynamics of proteins with an elastic network model. *Biophys J* 80:505–515.
- Bahar I, Atilgan AR, Erman B (1997) Direct evaluation of thermal fluctuations in proteins using a
  single-parameter harmonic potential. *Fold Des* 2:173–181.
- Bakan A, Bahar I (2009) The intrinsic dynamics of enzymes plays a dominant role in determining
  the structural changes induced upon inhibitor binding. *PNAS* 106:14349–14354.
- Bakan A, Meireles LM, Bahar I (2011) ProDy: protein dynamics inferred from theory and
  experiments. *Bioinformatics* 27:1575–1577.
- Brüschweiler R (1995) Collective protein dynamics and nuclear spin relaxation. *J Chem Phys*
  102:3396–3403.
- Gowers RJ *et al.* (2016) MDAnalysis: a Python package for the rapid analysis of molecular
  dynamics simulations. *Proc 15th Python in Science Conf*, 98–105.
- Kabsch W (1976) A solution for the best rotation to relate two sets of vectors. *Acta Cryst A*
  32:922–923.
- Kabsch W, Sander C (1983) Dictionary of protein secondary structure. *Biopolymers*
  22:2577–2637.
- Krüger DM, Ahmed A, Gohlke H (2012) NMSim web server: integrated approach for normal mode-based
  geometric simulations of biologically relevant conformational transitions in proteins.
  *Nucleic Acids Res* 40:W310–W316.
- Kurkcuoglu Z, Bahar I, Doruker P (2016) ClustENM: ENM-based sampling of essential conformational
  space at full atomic resolution. *J Chem Theory Comput* 12:4549–4562.
- Lezon TR, Bahar I (2010) Using entropy maximization to understand the determinants of structural
  dynamics beyond native contact topology. *PLoS Comput Biol* 6:e1000816.
- Marques O, Sanejouand YH (1995) Hinge-bending motion in citrate synthase arising from normal
  mode calculations. *Proteins* 23:557–560.
- Michaud-Agrawal N, Denning EJ, Woolf TB, Beckstein O (2011) MDAnalysis: a toolkit for the
  analysis of molecular dynamics simulations. *J Comput Chem* 32:2319–2327.
- Rego N, Koes D (2015) 3Dmol.js: molecular visualization with WebGL. *Bioinformatics*
  31:1322–1324.
- Tama F, Gadea FX, Marques O, Sanejouand YH (2000) Building-block approach for determining
  low-frequency normal modes of macromolecules. *Proteins* 41:1–7.
- Tama F, Sanejouand YH (2001) Conformational change of proteins arising from normal mode
  calculations. *Protein Eng* 14:1–6.
- Zhang S, Krieger JM, Zhang Y, Kaya C, Kaynak B, Mikulska-Ruminska K, Doruker P, Li H, Bahar I
  (2021) ProDy 2.0: increased scale and scope after 10 years of protein dynamics modelling with
  Python. *Bioinformatics* 37:3657–3659.
