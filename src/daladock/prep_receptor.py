"""
Stage 1 - Prepare the TARGET (receptor) -> receptor.pdbqt

prepare_receptor(inp, out) accepts:
  * a .pdb file
  * a 4-character PDB id (downloaded from RCSB)
  * a SMILES string (small-molecule 'receptor', built in 3D by RDKit)

Method is auto-selected: mostly standard amino acids -> Meeko, else OpenBabel.
Returns a dict of output paths + info.
"""
import os
import subprocess
import urllib.request

STD_AA = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
    "THR TRP TYR VAL MSE HID HIE HIP".split()
)
IONS = set("NA CL K MG ZN CA MN FE SO4 PO4 GOL EDO ACT".split())
VALID_TYPES = set("H HD HS C A N NA NS OA OS SA S P F Cl Br I Mg Zn Ca Mn Fe".split())


def _run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def _resolve_input(inp, base):
    """Return a path to a .pdb file for the receptor."""
    if os.path.isfile(inp):
        return inp
    if len(inp) == 4 and inp.isalnum():
        # use a local cache if present (HPC compute nodes often have no internet)
        for cache in (f"{inp.upper()}.pdb", os.path.join("..", "data", f"{inp.upper()}.pdb")):
            if os.path.isfile(cache):
                return cache
        url = f"https://files.rcsb.org/download/{inp.upper()}.pdb"
        dst = f"{inp.upper()}.pdb"
        urllib.request.urlretrieve(url, dst)
        return dst
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(inp)
    if m is None:
        raise ValueError(f"input '{inp}' is not a file, PDB id, or valid SMILES")
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(m)
    dst = f"{base}_fromsmiles.pdb"
    Chem.MolToPDBFile(m, dst)
    return dst


def _clean(src, base, keep_chains, strip_hetero, extract_ref):
    lines = open(src).read().splitlines()
    keep, hetres = [], {}
    for l in lines:
        rec = l[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        resn, ch = l[17:20].strip(), l[21]
        if resn == "HOH":
            continue
        if keep_chains and ch not in keep_chains:
            continue
        if rec == "HETATM":
            if resn in IONS:
                continue
            hetres.setdefault((resn, ch, l[22:26]), []).append(l)
            if strip_hetero:
                continue
        keep.append(l)
    clean = f"{base}_clean.pdb"
    open(clean, "w").write("\n".join(keep) + "\nEND\n")
    ref = None
    if extract_ref and hetres:
        biggest = max(hetres.values(), key=len)
        ref = f"{base}_ref.pdb"
        open(ref, "w").write("\n".join(biggest) + "\nEND\n")
    return clean, keep, ref


def _frac_standard(atomlines):
    # consider ALL kept residues (ATOM + HETATM); a protein is mostly standard AAs,
    # a glycopeptide like vancomycin is mostly non-standard -> picks the right prep tool
    res = {(l[17:20].strip(), l[21], l[22:26]) for l in atomlines
           if l[:6].strip() in ("ATOM", "HETATM")}
    return (sum(1 for r in res if r[0] in STD_AA) / len(res)) if res else 0.0


def _validate(pdbqt):
    types, zero, tot = set(), 0, 0
    for l in open(pdbqt):
        if l[:6].strip() in ("ATOM", "HETATM"):
            tot += 1
            types.add(l[77:79].strip())
            try:
                if abs(float(l[70:76])) < 1e-6:
                    zero += 1
            except ValueError:
                pass
    return {"types": sorted(types), "n_atoms": tot, "zero_charge": zero,
            "unknown_types": sorted(types - VALID_TYPES)}


def prepare_receptor(inp, out="receptor", method="auto", keep_chains=None,
                     strip_hetero=False, extract_ref=True, verbose=True):
    keep = set(keep_chains) if keep_chains else set()
    src = _resolve_input(inp, out)
    clean, atomlines, ref = _clean(src, out, keep, strip_hetero, extract_ref)

    if method == "auto":
        frac = _frac_standard(atomlines)
        method = "meeko" if frac > 0.5 else "obabel"
        if verbose:
            print(f"  standard-residue fraction={frac:.2f} -> method={method}")

    pdbqt = f"{out}.pdbqt"
    ok = False
    if method == "meeko":
        _run(f"mk_prepare_receptor.py --read_pdb {clean} -o {out} -p "
             f"--charge_model gasteiger -a")
        ok = os.path.isfile(pdbqt)
        if not ok and verbose:
            print("  Meeko failed; falling back to OpenBabel")
    if not ok:
        _run(f"obabel {clean} -O {pdbqt} -xr -h --partialcharge gasteiger")
        ok = os.path.isfile(pdbqt)
    if not ok:
        raise RuntimeError("receptor preparation failed")

    info = _validate(pdbqt)
    result = {"pdbqt": pdbqt, "clean_pdb": clean, "ref_ligand": ref,
              "method": method, **info}
    if verbose:
        print(f"  wrote {pdbqt}  types={' '.join(info['types'])}  "
              f"atoms={info['n_atoms']}")
        if ref:
            print(f"  extracted reference ligand -> {ref}")
        if info["unknown_types"]:
            print(f"  WARNING unknown atom types: {info['unknown_types']}")
    return result
