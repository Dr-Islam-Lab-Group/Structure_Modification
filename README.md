#RDKit Molecular Modification Pipeline
Overview
This repository contains a Python script that systematically generates chemical analogs from an input SDF file using SMARTS-based transformations and 3D structure generation with RDKit.
It is designed for medicinal chemistry optimization, bioisosteric replacement, and structure–property exploration.

The script performs two rounds of molecular transformations, generating optimized 3D conformations and saving all results as individual SDF files with a detailed transformation log.

#Key Features

==> SMARTS-driven functional group replacements

==> Two-round iterative molecular expansion

==> Bioisosteric, polarity, and lipophilicity-driven transformations

==> Automatic hydrogen addition, 3D embedding, and force-field optimization

==> Structured output directories for each round

==> CSV log of all transformations and generated molecules

#Transformation Categories

The script includes a curated library of medicinal chemistry transformations, including:

1.  Targeted group:
-CH3
    Modified groups:
-NH2, -NH(CH3), -N(CH3)2, -CF3, -CHO, -COOH, 
-CONH2, -CH2OH, -C6H4CN, -SCN, -SC6H5,
-C6H4NH2
2.  Targeted group:
-CN
    Modified groups:
-COOH, -CONH2, -NH2, -NH(CH3), -N(CH3)2, 
-NHC6H5, -C6H4CN, -CHO
3.  Targeted group:
    Gamma-lactam ring (2-pyrrolidinone)
4.  Modified groups
2-piperazinone
    Targeted group:
-CONH
    Modified groups
-NH2,-CONH2, -O�C(=O), �SO2?�NH2, -NH2, 
-NH(CH3), -N(CH3)2, -CHO, -COOH, -CONH2, -NHC6H5,-C6H4CN
5.  Targeted group:
-CF3
    Modified groups
-SO2NH2, -CN, -COOH, -CH3, -OCH3, -NH2, -CHO, -CONH2, 
-NH(CH3), -N(CH3)2, -NHC6H5, -C6H4CN
6.  Targeted group:
-F
    Modified groups
-NH2, -CHO, -COOH, CONH2, -NHCH3, 
-N(CH3)2, -C6H4NH2, -C6H4CN, -SCN, -SC6H5

#Requirements
1. Python Version

Python 3.8+

2. Required Libraries
pip install rdkit pandas


RDKit is recommended to be installed via conda:

conda install -c conda-forge rdkit

#Input

initial_molecules.sdf or Multi-molecule SDF file

#Output Files
1. SDF Files

==> One file per generated molecule

==> Fully hydrogenated

==> 3D embedded using ETKDG

==> Optimized with UFF

2. Transformation Log (transformation_log.csv)

Includes:

==> Transformation round

==> Parent molecule

==> Applied transformation

==> Product name

==> Output file path

#How It Works

==> Load molecules from the input SDF

==> Apply all SMARTS reactions (Round 1)

==> Generate 3D coordinates and optimize structures

==> Apply all reactions again to Round 1 products (Round 2)

==> Log all transformations to CSV

#Running the Script
python Structure_modification.py


Note: Ensure initial_molecules.sdf is in the same directory or update the path in the script.


If you use this workflow in research or publications, please cite RDKit:

RDKit: Open-source cheminformatics; http://www.rdkit.org

# Ligand Docking & MD Simulation Pipeline

This repository contains two stages of a structure-based drug discovery workflow:

1. **Docking** (`Docking_vina.sh`) — ligand preparation and virtual screening with AutoDock Vina.
2. **Molecular Dynamics** (`MD_inputs/`) — AMBER input files for minimization, equilibration, and production MD of a solvated system.

Documented in that order below, as requested. Note that the docking script's hardcoded Vina path (`.../After_MD_docking/Vina/vina`) suggests it may actually have been used *after* an MD run (e.g., docking against an MD-derived representative structure) rather than before — worth confirming which order applies to your project before running things.

---

## Requirements

| Stage | Software |
|---|---|
| Docking | [Open Babel](http://openbabel.org/) (`obabel`), [AutoDock Vina](http://vina.scripps.edu/) |
| MD | AmberTools / Amber (`sander` or `pmemd`, ideally `pmemd.cuda` for GPU) |

MD also requires a solvated, parameterized system (`prmtop` + `inpcrd`) built separately, e.g. with `tleap` or CHARMM-GUI — this is not included in `MD_inputs/`.

---

## 1. Docking (`Docking_vina.sh`)

### Expected directory layout

```
project/
├── conf.txt              # Vina config (you create this)
├── Docking_vina.sh
├── Ligands/
│   └── *.sdf              # input ligands (you provide these)
└── Results/                # created automatically
```

### Before running

- **Place ligand files** (`.sdf`) directly in `Ligands/`.
- **Create `conf.txt`** in the project root (one level above `Ligands/`). This is the Vina configuration file and must define the receptor and search box, e.g.:

  ```
  receptor = receptor.pdbqt
  center_x = 0.0
  center_y = 0.0
  center_z = 0.0
  size_x   = 20
  size_y   = 20
  size_z   = 20
  exhaustiveness = 8
  ```

  The receptor itself must be prepared as a `.pdbqt` beforehand (e.g. with AutoDockTools/MGLTools or `obabel`) — no receptor-prep step is included in the script.

- **Fix the hardcoded Vina path.** Line 31 currently points to:

  ```
  /Users/students-ad/Desktop/After_MD_docking/Vina/vina
  ```

  Change this to your local Vina binary path, or simply `vina` if it's on your `$PATH`.

### What the script does

1. Creates `Ligands/minimized/`, `Ligands/pdbqt/`, and `Results/`.
2. Energy-minimizes every `.sdf` in `Ligands/` with Open Babel (UFF force field, conjugate gradient, 200 steps, convergence criterion 0.1) → `Ligands/minimized/min_<name>.sdf`.
3. Converts each minimized ligand to PDBQT with Gasteiger partial charges (`obabel --gen3d --partialcharge gasteiger`) → `Ligands/pdbqt/<name>.pdbqt`.
4. Docks each ligand against the receptor with AutoDock Vina using `conf.txt` → poses to `Results/<name>.pdbqt`, full log to `Results/<name>.pdbqt.txt`.
5. Extracts the top-ranked (mode 1) binding affinity from each log and writes a summary to `results.csv` (`filename,affinity`).

### Running

```bash
chmod +x Docking_vina.sh
./Docking_vina.sh
```

### Output

- `Results/*.pdbqt` — docked poses
- `Results/*.pdbqt.txt` — full Vina logs (all binding modes and affinities)
- `results.csv` — summary table (filename, top-mode binding affinity in kcal/mol)

### Things to check before relying on this script

- The rename step `mv "$f" "${f%.pdbqt}.pdbqt"` currently renames each output file to itself — a no-op. If it was meant to strip some other suffix, that logic needs revisiting.
- `obabel` and `vina` must already be installed and reachable (aside from the hardcoded Vina binary path noted above).
- The script assumes one receptor (defined once in `conf.txt`) docked against many ligands — it doesn't support per-ligand receptor configs.

---

## 2. Molecular Dynamics (`MD_inputs/`)

Three AMBER (`sander`/`pmemd`) input decks implementing a standard minimize → NVT equilibrate → NPT production protocol in explicit solvent.

### Required (not included)

- `prmtop` (topology) and `inpcrd` (coordinates) for the solvated system — build with `tleap`/AmberTools or CHARMM-GUI.
- The same starting coordinates are used as the reference structure for positional restraints (`ntr=1`) in all three steps.

### Protocol

| File | Ensemble | Purpose | Length | Timestep |
|---|---|---|---|---|
| `step4.0_minimization.mdin` | — | Energy minimization, restraints on (5000 steps: 2500 steepest-descent → conjugate gradient) | n/a | n/a |
| `step4.1_equilibration.mdin` | NVT | Equilibration, Langevin thermostat at 303.15 K, restraints on | 2,000,000 steps = 2 ns | 1 fs |
| `step5_production.mdin` | NPT | Production, Langevin thermostat (303.15 K) + Monte Carlo barostat (1 bar), restraints on | 125,000,000 steps = 250 ns | 2 fs |

Shared settings: 12 Å nonbonded cutoff with force-based switching from 10 Å, SHAKE on hydrogen-containing bonds (equilibration/production), NetCDF trajectory and restart output, water residues named `WAT`.

**Positional restraints:** all three steps use `ntr=1` with mask `RES 1 299 300 300` at 1.0 kcal/mol·Å². This restrains a specific residue range (protein/ligand) — **you must edit this to match your own system's residue numbering** before running; it will not transfer as-is to a different prmtop.

### Running (example with `pmemd.cuda`; substitute `sander` for CPU-only runs)

```bash
# Minimization
pmemd.cuda -O -i step4.0_minimization.mdin -o min.out \
  -p system.prmtop -c system.inpcrd -r min.rst7 -ref system.inpcrd

# NVT equilibration
pmemd.cuda -O -i step4.1_equilibration.mdin -o eq.out \
  -p system.prmtop -c min.rst7 -r eq.rst7 -x eq.nc -ref min.rst7

# NPT production
pmemd.cuda -O -i step5_production.mdin -o prod.out \
  -p system.prmtop -c eq.rst7 -r prod.rst7 -x prod.nc -ref eq.rst7
```

### Output

- `*.out` — energy/log files
- `*.nc` — NetCDF trajectories
- `*.rst7` — restart/coordinate files for chaining to the next step

### Things to check before relying on these inputs

- Production (`step5`) is configured for **250 ns** in one run (`nstlim=125000000`, `dt=0.002`) — confirm this matches your intended run length, and plan for checkpointing/restarts given the wall-clock this implies.
- The restraint mask `RES 1 299 300 300` is system-specific and must be verified against your `prmtop` before use.
- No `step6`-style graded restraint release is included; restraints stay on (`ntr=1`, full weight) through equilibration and production as written.

---


