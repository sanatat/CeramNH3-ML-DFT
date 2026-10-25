from pathlib import Path
from ase.io import read

base = Path("outputs")
source = base / "qe_Cu3Fe9_clean"
structures = base / "nrr_structures/Cu3Fe9"

n_pseudo = Path.home() / "qe-source/qe-7.5/pseudo/N-PBE.upf"

for site in ["Fe", "Cu"]:
    out = base / f"qe_Cu3Fe9_N2_{site}"
    pseudo = out / "pseudo"
    pseudo.mkdir(parents=True, exist_ok=True)
    (out / "tmp").mkdir(exist_ok=True)

    for file in (source / "pseudo").glob("*.UPF"):
        (pseudo / file.name).write_bytes(file.read_bytes())
    (pseudo / n_pseudo.name).write_bytes(n_pseudo.read_bytes())

    atoms = read(structures / f"N2_on_{site}_unrelaxed.traj")

    text = [
        "&CONTROL",
        " calculation = 'scf',",
        f" prefix = 'Cu3Fe9_N2_{site}',",
        " outdir = './tmp',",
        " pseudo_dir = './pseudo',",
        "/",
        "&SYSTEM",
        " ibrav = 0,",
        f" nat = {len(atoms)},",
        " ntyp = 3,",
        " ecutwfc = 70,",
        " ecutrho = 560,",
        " nspin = 2,",
        " starting_magnetization(1) = 0.35,",
        " starting_magnetization(2) = 0.02,",
        " starting_magnetization(3) = 0.0,",
        " occupations = 'smearing',",
        " smearing = 'mv',",
        " degauss = 0.02,",
        "/",
        "&ELECTRONS",
        " conv_thr = 1.0d-6,",
        " electron_maxstep = 150,",
        " mixing_beta = 0.2,",
        "/",
        "ATOMIC_SPECIES",
        "Fe 55.845 Fe.jry.pbe.UPF",
        "Cu 63.546 Cu.pbe-dn-kjpaw_psl.0.2.UPF",
        "N 14.007 N-PBE.upf",
        "CELL_PARAMETERS angstrom"
    ]

    for v in atoms.cell.array:
        text.append(" ".join(f"{x:.10f}" for x in v))

    text.append("ATOMIC_POSITIONS angstrom")
    for atom in atoms:
        x, y, z = atom.position
        text.append(f"{atom.symbol} {x:.10f} {y:.10f} {z:.10f}")

    text += ["K_POINTS automatic", "4 4 1 0 0 0"]

    path = out / "n2_scf.in"
    path.write_text("\n".join(text) + "\n")
    print("Created:", path)
