# dala-dock

An interactive teaching pipeline for **molecular docking**, built around a real story:
**why the antibiotic vancomycin stops working against resistant bacteria.**

You dock small peptide "candidates" into vancomycin and watch, stage by stage, why the
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
   git clone https://github.com/<your-username>/dala-dock.git
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
git clone https://github.com/<your-username>/dala-dock.git
cd dala-dock
conda env create -f environment.yml        # laptop env (Vina, CPU)
conda activate dala-dock
jupyter lab notebooks/docking_demo.ipynb
```
Same notebook; cell 1 prints `engine ... vina`. Scores differ from AutoDock-GPU (different
scoring functions) — compare rankings *within* one engine, not absolute numbers across engines.

---

## For the maintainer: build the shared environment on Perlmutter

Done **once** by the project owner; the coordinator then makes it available to participants.

1. Clone the repo on Perlmutter and edit the project number in
   `scripts/build_on_perlmutter.sh` and `kernels/dala-dock/kernel-helper.sh`
   (`$CFS/m4833/projects/project<N>`).
2. Build the shared env + AutoDock-GPU into CFS:
   ```bash
   bash scripts/build_on_perlmutter.sh
   ```
   This creates `$CFS/m4833/projects/project<N>/env` and builds AutoDock-GPU for the A100.
   It prints the `PROJECT_ROOT` and `DALADOCK_CUDA_LIB` values for the kernel.
3. Put those values into `kernels/dala-dock/kernel-helper.sh`, then register the kernel:
   ```bash
   mkdir -p $HOME/.local/share/jupyter/kernels/dala-dock
   cp kernels/dala-dock/kernel.json kernels/dala-dock/kernel-helper.sh \
      $HOME/.local/share/jupyter/kernels/dala-dock/
   chmod +x $HOME/.local/share/jupyter/kernels/dala-dock/kernel-helper.sh
   ```
4. Test `notebooks/docking_demo.ipynb` on jupyter.nersc.gov (Perlmutter GPU) with the
   `dala-dock` kernel, then share the env path (`$CFS/m4833/projects/project<N>/env`)
   with the coordinator.

See `configs/nersc.yaml` for the account/partition/module settings.

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

## Caveats to teach
- Empirical scoring → **trends are more trustworthy than absolute kcal/mol**.
- Docking predicts **binding**, not antibiotic efficacy (binding is necessary, not sufficient).
- Vina and AutoDock-GPU use different scoring functions — compare rankings within one engine.
- Automatic pocket detection can pick the wrong cavity — always sanity-check the box.
