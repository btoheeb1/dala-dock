"""
Visualization helpers (py3Dmol) so each stage can be *seen* in the notebook.

  view_structure(path)                 - show a receptor/molecule
  add_box(view, center, size)          - draw the search box
  view_complex(receptor, pose)         - receptor + docked ligand pose
  view_ligand(path)                    - a single ligand as sticks
"""
import subprocess
import py3Dmol


def _to_pdb(path, first_only=False):
    """Convert any OpenBabel-readable file to a PDB string (for 3Dmol)."""
    cmd = f"obabel {path} -opdb" + (" -f 1 -l 1" if first_only else "")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    return out


def _has_protein(pdb_text):
    aa = {"ALA", "GLY", "SER", "LEU", "VAL", "THR", "LYS", "ASP", "GLU", "PHE"}
    return any(l[17:20].strip() in aa for l in pdb_text.splitlines()
              if l[:6].strip() in ("ATOM", "HETATM"))


def view_structure(path, width=640, height=480, style="auto", surface=False):
    pdb = _to_pdb(path)
    v = py3Dmol.view(width=width, height=height)
    v.addModel(pdb, "pdb")
    if style == "auto":
        style = "cartoon" if _has_protein(pdb) else "stick"
    if style == "cartoon":
        v.setStyle({"cartoon": {"color": "spectrum"}})
        v.addStyle({"hetflag": True}, {"stick": {"colorscheme": "greenCarbon"}})
    else:
        v.setStyle({"stick": {"colorscheme": "cyanCarbon"}})
    if surface:
        v.addSurface(py3Dmol.VDW, {"opacity": 0.6, "color": "white"})
    v.zoomTo()
    return v


def add_box(view, center, size, color="magenta", opacity=0.35, wireframe=True):
    view.addBox({
        "center": {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])},
        "dimensions": {"w": float(size[0]), "h": float(size[1]), "d": float(size[2])},
        "color": color, "opacity": opacity, "wireframe": wireframe,
    })
    return view


def view_box(receptor, center, size, width=640, height=480, style="auto"):
    v = view_structure(receptor, width, height, style)
    add_box(v, center, size)
    return v


def view_ligand(path, width=480, height=360):
    v = py3Dmol.view(width=width, height=height)
    v.addModel(_to_pdb(path, first_only=True), "pdb")
    v.setStyle({"stick": {"colorscheme": "cyanCarbon"}})
    v.zoomTo()
    return v


def _split_models(poses_pdbqt):
    """Return a list of single-model PDB strings from a multi-pose file."""
    pdb = _to_pdb(poses_pdbqt)
    blocks = pdb.split("ENDMDL")
    return [b + "ENDMDL\n" for b in blocks if ("ATOM" in b or "HETATM" in b)]


_POSE_COLORS = ["red", "orange", "yellow", "green", "cyan", "blue",
                "purple", "magenta", "lime", "white"]


def view_poses(receptor, poses_pdbqt, max_poses=10, width=780, height=560,
               receptor_style="auto"):
    """Overlay all docked poses (each a different color) inside the receptor."""
    rec_pdb = _to_pdb(receptor)
    v = py3Dmol.view(width=width, height=height)
    v.addModel(rec_pdb, "pdb")
    if receptor_style == "auto":
        receptor_style = "cartoon" if _has_protein(rec_pdb) else "stick"
    if receptor_style == "cartoon":
        v.setStyle({"cartoon": {"color": "lightgray"}})
    else:
        v.setStyle({"stick": {"colorscheme": "grayCarbon", "radius": 0.08}})
    for i, m in enumerate(_split_models(poses_pdbqt)[:max_poses]):
        v.addModel(m, "pdb")
        v.setStyle({"model": -1},
                   {"stick": {"color": _POSE_COLORS[i % len(_POSE_COLORS)], "radius": 0.12}})
    v.zoomTo({"model": -1})
    return v


def view_pose_n(receptor, poses_pdbqt, n=1, label=None, width=720, height=520,
                receptor_style="auto"):
    """Show a single docked pose (1-indexed) inside the receptor."""
    rec_pdb = _to_pdb(receptor)
    models = _split_models(poses_pdbqt)
    idx = max(1, min(int(n), len(models))) - 1
    v = py3Dmol.view(width=width, height=height)
    v.addModel(rec_pdb, "pdb")
    if receptor_style == "auto":
        receptor_style = "cartoon" if _has_protein(rec_pdb) else "stick"
    if receptor_style == "cartoon":
        v.setStyle({"cartoon": {"color": "lightgray"}})
    else:
        v.setStyle({"stick": {"colorscheme": "grayCarbon", "radius": 0.08}})
    v.addModel(models[idx], "pdb")
    v.setStyle({"model": -1}, {"stick": {"colorscheme": "orangeCarbon", "radius": 0.22}})
    if label:
        v.addLabel(str(label), {"backgroundColor": "black", "fontColor": "white",
                                "fontSize": 12})
    v.zoomTo({"model": -1})
    return v


def _no_atoms(lines):
    """Electronegative (N/O) heavy atoms -> list of (x,y,z); H-bond donors/acceptors."""
    out = []
    for l in lines:
        if l[:6].strip() in ("ATOM", "HETATM"):
            t = l[77:79].strip()
            if t[:1] in ("N", "O"):
                out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def _model_lines(poses_pdbqt, n):
    """ATOM/HETATM lines of the n-th pose (1-indexed) from a multi-model file."""
    models, cur = [], []
    for l in open(poses_pdbqt):
        if l.startswith("MODEL"):
            cur = []
        elif l.startswith("ENDMDL"):
            if cur:
                models.append(cur); cur = []
        elif l[:6].strip() in ("ATOM", "HETATM"):
            cur.append(l)
    if cur:
        models.append(cur)
    if not models:
        models = [[l for l in open(poses_pdbqt) if l[:6].strip() in ("ATOM", "HETATM")]]
    return models[max(1, min(int(n), len(models))) - 1]


def view_hbonds(receptor, poses_pdbqt, pose=1, cutoff=3.5, min_dist=2.4,
                width=800, height=600, receptor_style="stick"):
    """Show a docked pose with candidate<->receptor hydrogen bonds as dashed lines.

    H-bonds are estimated the common quick way: pairs of electronegative atoms
    (N/O on the ligand and N/O on the receptor) whose distance falls in the
    H-bonding range (~2.4-3.5 A). These dashed lines are the *non-covalent*
    interactions that hold the candidate in the pocket - not chemical bonds.
    """
    import math
    rec_lines = [l for l in open(receptor) if l[:6].strip() in ("ATOM", "HETATM")]
    lig_lines = _model_lines(poses_pdbqt, pose)
    rec_no, lig_no = _no_atoms(rec_lines), _no_atoms(lig_lines)

    v = py3Dmol.view(width=width, height=height)
    v.addModel(_to_pdb(receptor), "pdb")
    if receptor_style == "surface":
        v.setStyle({"stick": {"colorscheme": "grayCarbon", "radius": 0.08}})
        v.addSurface(py3Dmol.VDW, {"opacity": 0.4, "color": "white"})
    else:
        v.setStyle({"stick": {"colorscheme": "grayCarbon", "radius": 0.1}})
    v.addModel(_split_models(poses_pdbqt)[max(1, int(pose)) - 1], "pdb")
    v.setStyle({"model": -1}, {"stick": {"colorscheme": "orangeCarbon", "radius": 0.2}})

    n = 0
    for lx, ly, lz in lig_no:
        for rx, ry, rz in rec_no:
            d = math.dist((lx, ly, lz), (rx, ry, rz))
            if min_dist <= d <= cutoff:
                v.addCylinder({"start": {"x": lx, "y": ly, "z": lz},
                               "end": {"x": rx, "y": ry, "z": rz},
                               "radius": 0.06, "color": "yellow",
                               "dashed": True, "fromCap": 1, "toCap": 1})
                v.addLabel(f"{d:.1f}", {"position": {"x": (lx + rx) / 2,
                                        "y": (ly + ry) / 2, "z": (lz + rz) / 2},
                                        "fontSize": 10, "fontColor": "black",
                                        "backgroundColor": "white", "backgroundOpacity": 0.6})
                n += 1
    print(f"  {n} candidate<->receptor H-bond contacts (<= {cutoff} A) shown as dashed lines")
    v.zoomTo({"model": -1})
    return v


def view_pose_vs_reference(receptor, poses_pdbqt, ref_pdb, pose=1,
                           width=760, height=560, receptor_style="auto"):
    """Overlay a docked pose (orange) with a known reference ligand (green) in the
    receptor - a visual validation that docking landed in the right place."""
    rec_pdb = _to_pdb(receptor)
    v = py3Dmol.view(width=width, height=height)
    v.addModel(rec_pdb, "pdb")
    if receptor_style == "auto":
        receptor_style = "cartoon" if _has_protein(rec_pdb) else "stick"
    if receptor_style == "cartoon":
        v.setStyle({"cartoon": {"color": "lightgray"}})
    else:
        v.setStyle({"stick": {"colorscheme": "grayCarbon", "radius": 0.08}})
    v.addModel(_split_models(poses_pdbqt)[max(1, int(pose)) - 1], "pdb")   # docked pose
    v.setStyle({"model": -1}, {"stick": {"colorscheme": "orangeCarbon", "radius": 0.2}})
    v.addModel(_to_pdb(ref_pdb), "pdb")                                    # reference
    v.setStyle({"model": -1}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}})
    v.zoomTo({"model": -1})
    return v


def view_complex(receptor, pose, width=700, height=520, receptor_style="auto",
                 surface=False, pose_color="orangeCarbon"):
    rec_pdb = _to_pdb(receptor)
    lig_pdb = _to_pdb(pose, first_only=True)   # best (first) pose
    v = py3Dmol.view(width=width, height=height)
    v.addModel(rec_pdb, "pdb")
    if receptor_style == "auto":
        receptor_style = "cartoon" if _has_protein(rec_pdb) else "stick"
    if receptor_style == "cartoon":
        v.setStyle({"cartoon": {"color": "lightgray"}})
        v.addStyle({"hetflag": True}, {"stick": {"colorscheme": "grayCarbon"}})
    else:
        v.setStyle({"stick": {"colorscheme": "grayCarbon", "radius": 0.12}})
    if surface:
        v.addSurface(py3Dmol.VDW, {"opacity": 0.5, "color": "white"})
    v.addModel(lig_pdb, "pdb")
    v.setStyle({"model": -1}, {"stick": {"colorscheme": pose_color, "radius": 0.2}})
    v.zoomTo({"model": -1})
    return v
