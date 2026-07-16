"""
Stage 2 - Define the search BOX (where docking is allowed to look).

define_box(receptor_pdbqt, ...) returns dict {center:[x,y,z], size:[sx,sy,sz]}.

Detection modes:
  ref     : center + size from a reference ligand PDB (most reliable)
  fpocket : run fpocket, use the top-ranked pocket (verify it! can be wrong)
  blind   : box encloses the whole receptor
  auto    : ref if given, else fpocket, else blind
Manual override: pass center=[...] and/or size=[...].
"""
import glob
import os
import subprocess
import numpy as np


def _run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def _coords(path):
    P = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
         for l in open(path) if l[:6].strip() in ("ATOM", "HETATM")]
    return np.array(P)


def _box_from_points(P, pad):
    return P.mean(0), (P.max(0) - P.min(0)) + pad


def _fpocket_box(receptor_pdbqt, base, pad):
    pdb = f"{base}_forfpocket.pdb"
    _run(f"obabel {receptor_pdbqt} -O {pdb}")
    _run(f"fpocket -f {pdb}")
    vert = sorted(glob.glob(f"{pdb[:-4]}_out/pockets/pocket*_vert.pqr"))
    if not vert:
        return None
    P = _coords(vert[0])
    return _box_from_points(P, pad) if len(P) else None


def define_box(receptor_pdbqt, ref=None, detect="auto", center=None, size=None,
               pad=8.0, base=None, verbose=True):
    base = base or os.path.splitext(receptor_pdbqt)[0]
    c = np.array(center, float) if center is not None else None
    s = np.array(size, float) if size is not None else None

    if c is None:
        order = {"auto": ["ref", "fpocket", "blind"], "ref": ["ref"],
                 "fpocket": ["fpocket"], "blind": ["blind"]}[detect]
        for mode in order:
            if mode == "ref" and ref and os.path.isfile(ref):
                c, s2 = _box_from_points(_coords(ref), pad)
                s = s if s is not None else s2
                src = f"reference ligand {ref}"
                break
            if mode == "fpocket":
                res = _fpocket_box(receptor_pdbqt, base, pad)
                if res:
                    c, s2 = res
                    s = s if s is not None else s2
                    src = "fpocket (VERIFY location!)"
                    break
            if mode == "blind":
                c, s2 = _box_from_points(_coords(receptor_pdbqt), pad)
                s = s if s is not None else s2
                src = "blind (whole receptor)"
                break
    else:
        src = "manual"
        if s is None:
            s = np.array([22.5, 22.5, 22.5])

    if c is None:
        raise RuntimeError("could not determine a box")
    result = {"center": [round(float(x), 3) for x in c],
              "size": [round(float(x), 1) for x in s], "source": src}
    if verbose:
        print(f"  box from {src}: center={result['center']} size={result['size']}")
    return result
