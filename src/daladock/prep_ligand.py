"""
Stage 3 - Prepare a CANDIDATE (ligand) -> <out>.sdf + <out>.pdbqt

prepare_ligand(inp, out) accepts a SMILES string or a file (.sdf/.mol/.mol2/.pdb).
Steps: RDKit (add H, embed 3D, MMFF) -> SDF -> Meeko -> PDBQT.
"""
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem


def _run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def _load(inp):
    if os.path.isfile(inp):
        ext = os.path.splitext(inp)[1].lower()
        loader = {".sdf": lambda p: next(Chem.SDMolSupplier(p, removeHs=False)),
                  ".mol": lambda p: Chem.MolFromMolFile(p, removeHs=False),
                  ".mol2": lambda p: Chem.MolFromMol2File(p, removeHs=False),
                  ".pdb": lambda p: Chem.MolFromPDBFile(p, removeHs=False)}.get(ext)
        if loader is None:
            raise ValueError(f"unsupported file type {ext}")
        m = loader(inp)
        has3d = m is not None and m.GetNumConformers() > 0
    else:
        m, has3d = Chem.MolFromSmiles(inp), False
    if m is None:
        raise ValueError(f"could not parse ligand '{inp}'")
    return m, has3d


def prepare_ligand(inp, out, seed=42, verbose=True):
    m, has3d = _load(inp)
    m = Chem.AddHs(m, addCoords=has3d)
    if not has3d:
        AllChem.EmbedMolecule(m, randomSeed=seed)
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    sdf, pdbqt = f"{out}.sdf", f"{out}.pdbqt"
    Chem.MolToMolFile(m, sdf)
    # Use Meeko's Python API in-process (avoids PATH/subprocess picking a wrong meeko)
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    setups = MoleculePreparation().prepare(m)
    if not setups:
        raise RuntimeError("Meeko could not prepare this ligand")
    written = PDBQTWriterLegacy.write_string(setups[0])
    pdbqt_string = written[0] if isinstance(written, tuple) else written
    if not pdbqt_string or not str(pdbqt_string).strip():
        raise RuntimeError(f"Meeko produced empty PDBQT (write result: {written!r})")
    with open(pdbqt, "w") as fh:
        fh.write(pdbqt_string)
    types = sorted({l[77:79].strip() for l in open(pdbqt)
                    if l[:6].strip() in ("ATOM", "HETATM")})
    tors = [l for l in open(pdbqt) if l.startswith("REMARK") and "active torsions" in l]
    ntors = int(tors[0].split()[1]) if tors else 0
    if verbose:
        print(f"  {out}: {m.GetNumAtoms()} atoms, {ntors} torsions, types={' '.join(types)}")
    return {"sdf": sdf, "pdbqt": pdbqt, "types": types, "n_torsions": ntors,
            "smiles": Chem.MolToSmiles(Chem.RemoveHs(m))}
