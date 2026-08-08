# dala-dock

An interactive teaching pipeline for **molecular docking**.

You dock small peptide "candidates" (e.g **D-Ala-D-Ala**) into a target (e.g vancomycin) and watch, stage by stage, why the
native **D-Ala-D-Ala** cell-wall terminus binds tightly while the resistant
**D-Ala-D-Lac** variant does not.

Same notebook runs two ways — the engine is chosen automatically (`engine="auto"`):
- **GPU node (NERSC Perlmutter)** → **AutoDock-GPU** on an A100.
- **Laptop / no GPU** → **AutoDock Vina** (CPU). Nothing else needed.

---

## For students (NERSC Perlmutter)

You do **not** install anything — a shared environment is provided.

1. Log in to **https://jupyter.nersc.gov** and start a server on a **Perlmutter GPU** node.
2. Open a terminal in Jupyter and clone the repo into your space:
   ```bash
   cd $SCRATCH        # or your preferred directory
   git clone https://github.com/btoheeb1/dala-dock.git
   cd dala-dock
   ```
3. Open `notebooks/docking_demo.ipynb`.
4. Top-right, select the **`dala-dock`** kernel (the shared environment).
5. **Run All** (or step through cell by cell). Cell 1 should print `engine ... adgpu`.

Everything (RDKit, Meeko, Vina, autogrid, AutoDock-GPU, py3Dmol) is already in the shared
kernel — you just run the notebook.

---

## What the notebook walks through

| Stage | Module | What you see |
|---|---|---|
| 1. Prepare receptor | `daladock.prep_receptor` | vancomycin from the PDB, in 3D |
| 2. Define the box | `daladock.box` | the search box drawn over the binding site |
| 3. Prepare candidates | `daladock.prep_ligand` | each D-Ala peptide built from SMILES |
| 4. Dock | `daladock.dock` | each candidate docked (Vina or AutoDock-GPU) |
| 4b. Explore poses | `daladock.analyze` + `viz` | per-pose table, clustering, all poses in 3D, best vs worst |
| 5. Analyze | `daladock.analyze` | ranked leaderboard (sensitive vs resistant) |

Edit `data/candidates.csv` to try your own molecules.

---

## Run it on your own laptop (no GPU, uses Vina)

```bash
git clone https://github.com/btoheeb1/dala-dock.git
cd dala-dock
conda env create -f environment.yml        # laptop env (Vina, CPU)
conda activate dala-dock
jupyter lab notebooks/docking_demo.ipynb
```
Same notebook; cell 1 prints `engine ... vina`. Scores differ from AutoDock-GPU (different
scoring functions) — compare rankings *within* one engine, not absolute numbers across engines.

## Exercises (student homework)

Do these in the **"Your turn"** section at the end of `notebooks/docking_demo.ipynb`.

**Tier 1 — Interpret + validate**
- From the leaderboard, which terminus binds tightest? Does the order match the resistance
  story (native **D-Ala-D-Ala** tighter than resistant **D-Ala-D-Lac** / **D-Ala-D-Ser**)?
- In Stage 4b, how many poses fall in the *best* cluster? What does a large low-energy
  cluster tell you about confidence?
- **Validate:** overlay the best docked pose with the crystal reference (`data/dala_ref.pdb`)
  and check it lands in the real pocket. Then argue: why is *binding* necessary but not
  sufficient for an antibiotic?

**Tier 2 — Tweak the parameters**
- Increase the number of poses (`n_poses` for Vina, `nrun` for AutoDock-GPU) and
  `exhaustiveness`. Does the top pose or its cluster change? Note the runtime — that's the
  **quality-vs-cost trade-off** at the heart of HPC.

**Tier 3 — Add & screen candidates**
- Add your own molecule(s) to `data/candidates.csv` (find real SMILES on
  [PubChem](https://pubchem.ncbi.nlm.nih.gov) and record the CID). Predict, then test.
- **Scale up:** screen the ~50 PubChem-verified molecules in `data/candidates_scan.csv`
  and build the leaderboard. Which terminal residues does vancomycin tolerate?

---

## Using your own receptor / candidates

The pipeline is general:
1. **Receptor:** give `prepare_receptor()` a PDB file, a 4-letter PDB id, or a SMILES.
2. **Box:** center on a bound ligand (`detect="ref"`), let `fpocket` try (`detect="fpocket"`
   — *verify it*), or search the whole target (`detect="blind"`).
3. **Candidates:** put SMILES in a CSV and loop `prepare_ligand()` + `dock()`.

**Validate first:** if your target has a known bound ligand, redock it and confirm the pose
lands within ~2 Å of the crystal before trusting scores on new molecules.

---

## Layout
```
environment.yml            laptop env (Vina, CPU)
environment_nersc.yml      shared HPC env (adds autogrid, fpocket, ipykernel) for CFS
configs/nersc.yaml         Perlmutter account / partition / modules
configs/loni.yaml          example config for another cluster
kernels/dala-dock/         Jupyter kernel (kernel.json + helper) for the shared env
scripts/build_on_perlmutter.sh   one-time build of the shared env + AutoDock-GPU
data/candidates.csv        the D-Ala candidate panel (edit me)
data/1FVM.pdb              cached structure (works with no internet on compute nodes)
src/daladock/              the pipeline (one module per stage) + viz helpers
notebooks/docking_demo.ipynb   the interactive tutorial
```

## Data & attribution

- **Receptor + reference ligand:** vancomycin and the bound D-Ala-D-Ala come from
  **RCSB PDB entry [1FVM](https://www.rcsb.org/structure/1FVM)** (1.8 Å complex of vancomycin
  with di-acetyl-Lys-D-Ala-D-Ala), downloaded from
  https://files.rcsb.org/download/1FVM.pdb (cached in `data/1FVM.pdb`).
  Citation: Nitanai Y, Kikuchi T, Kakoi K, Hanamaki S, Fujisawa I, Aoki K.
  *"Crystal Structures of the Complexes between Vancomycin and Cell-Wall Precursor Analogs."*
  **J. Mol. Biol.** 385(5):1422–1442 (2009). doi:10.1016/j.jmb.2008.10.026 (PMID 18976660).
- **Candidate molecules:** sourced from **PubChem** (https://pubchem.ncbi.nlm.nih.gov).
  Every row in `data/candidates.csv` and `data/candidates_scan.csv` records the molecule's
  **PubChem CID, InChIKey, and source URL** — nothing is hand-fabricated. 3D coordinates are
  generated from those SMILES with RDKit at run time.
- RCSB PDB and PubChem data are freely available for research and education.


