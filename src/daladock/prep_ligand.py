"""
Stage 3 - Prepare a CANDIDATE (ligand) -> <out>.sdf + <out>.pdbqt

prepare_ligand(inp, out) accepts a SMILES string or a file (.sdf/.mol/.mol2/.pdb).
Steps: RDKit (add H, embed 3D, MMFF) -> SDF -> Meeko -> PDBQT.
"""
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

STD_AA = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
    "THR TRP TYR VAL MSE HID HIE HIP".split()
)


def _run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def _pdb_is_smallmol(path):
    """A peptide/protein PDB has standard-AA residues (RDKit knows their bonds);
    a small-molecule ligand PDB does not, and carries no bond-order info."""
    res = {(l[17:20].strip(), l[21], l[22:26]) for l in open(path)
           if l[:6].strip() in ("ATOM", "HETATM")}
    return sum(1 for r in res if r[0] in STD_AA) < 2


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


def _prep_smallmol_pdb_obabel(inp, out, rigid, verbose):
    """A small-molecule ligand PDB carries no bond-order info, so RDKit can't
    rebuild it (mis-perceives bonds -> rigid, or fails on charges). OpenBabel
    perceives bonds from geometry and writes a valid PDBQT directly."""
    sdf, pdbqt = f"{out}.sdf", f"{out}.pdbqt"
    _run(f"obabel {inp} -O {sdf}")                       # keep an SDF for inspection
    flags = "-h --partialcharge gasteiger" + (" -xr" if rigid else "")
    r = _run(f"obabel {inp} -O {pdbqt} {flags}")
    if not os.path.isfile(pdbqt) or not open(pdbqt).read().strip():
        raise RuntimeError(f"OpenBabel could not prepare ligand '{inp}':\n{r.stderr[-300:]}")
    atoms = [l for l in open(pdbqt) if l[:6].strip() in ("ATOM", "HETATM")]
    types = sorted({l[77:79].strip() for l in atoms})
    tdof = [l for l in open(pdbqt) if l.startswith("TORSDOF")]
    ntors = 0 if rigid else (int(tdof[0].split()[1]) if tdof else
                             sum(1 for l in open(pdbqt) if l.startswith("BRANCH")))
    smi = _run(f"obabel {inp} -osmi 2>/dev/null").stdout.split()
    if verbose:
        print(f"  {out}: {len(atoms)} atoms, {ntors} torsions, "
              f"types={' '.join(types)}  (OpenBabel)")
    return {"sdf": sdf, "pdbqt": pdbqt, "types": types, "n_torsions": ntors,
            "smiles": smi[0] if smi else None}


def _flatten_rigid(pdbqt_string):
    """Collapse the torsion tree into one rigid body: all atoms in a single ROOT,
    no BRANCHes, TORSDOF 0. For large molecules (e.g. peptides) that have far more
    rotatable bonds than a docking engine can handle."""
    atoms, remarks = [], []
    for l in pdbqt_string.splitlines():
        if l[:6].strip() in ("ATOM", "HETATM"):
            atoms.append(l)
        elif l.startswith("REMARK") and "torsion" not in l.lower():
            remarks.append(l)
    return "\n".join(remarks + ["ROOT"] + atoms + ["ENDROOT", "TORSDOF 0"]) + "\n"


def prepare_ligand(inp, out, seed=42, rigid=False, verbose=True):
    # Small-molecule ligand given as a bare PDB: RDKit can't rebuild bonds from it,
    # so route through OpenBabel (geometry-based bond perception) straight to PDBQT.
    if (os.path.isfile(inp) and os.path.splitext(inp)[1].lower() == ".pdb"
            and _pdb_is_smallmol(inp)):
        return _prep_smallmol_pdb_obabel(inp, out, rigid, verbose)
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
    # rigid=True: keep disulfide/macrocyclic rings intact (no CG/G glue atoms that
    # autogrid cannot map) and dock the whole molecule as a single rigid body.
    setups = MoleculePreparation(rigid_macrocycles=rigid).prepare(m)
    if not setups:
        raise RuntimeError("Meeko could not prepare this ligand")
    written = PDBQTWriterLegacy.write_string(setups[0])
    pdbqt_string = written[0] if isinstance(written, tuple) else written
    if not pdbqt_string or not str(pdbqt_string).strip():
        raise RuntimeError(f"Meeko produced empty PDBQT (write result: {written!r})")
    if rigid:
        pdbqt_string = _flatten_rigid(pdbqt_string)
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
