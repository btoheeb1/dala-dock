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
