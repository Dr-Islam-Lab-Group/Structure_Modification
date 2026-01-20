from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import os
import pandas as pd

input_sdf = "initial_molecules.sdf"  # Replace with your actual file name
round1_dir = "round1_outputs_3D"
round2_dir = "round2_outputs_3D"
os.makedirs(round1_dir, exist_ok=True)
os.makedirs(round2_dir, exist_ok=True)

log_entries = []  # Store logs for CSV

# SMARTS DICTIONARY
# SMARTS LIBRARY FROM YOUR TABLE
reaction_smarts_dict = {

    # 1) Target: -CH3
    "methyl → NH2":          "[CH3:1]>>[CH2:1]N",
    "methyl → NHCH3":        "[CH3:1]>>[CH2:1]N(C)",
    "methyl → N(CH3)2":      "[CH3:1]>>[CH2:1]N(C)C",
    "methyl → CF3":          "[CH3:1]>>[C:1](F)(F)F",
    "methyl → CHO":          "[CH3:1]>>[CH:1]=O",
    "methyl → COOH":         "[CH3:1]>>[C:1](=O)O",
    "methyl → CONH2":        "[CH3:1]>>[C:1](=O)N",
    "methyl → CH2OH":        "[CH3:1]>>[CH2:1]O",
    "methyl → C6H4CN":       "[CH3:1]>>[CH2:1]c1ccc(C#N)cc1",
    "methyl → SCN":          "[CH3:1]>>[CH2:1]SC#N",
    "methyl → SC6H5":        "[CH3:1]>>[CH2:1]Sc1ccccc1",
    "methyl → C6H4NH2":      "[CH3:1]>>[CH2:1]c1ccc(N)cc1",

    # 2) Target: -CN
    "nitrile → COOH":        "[C:1]#[N]>>[C:1](=O)O",
    "nitrile → CONH2":       "[C:1]#[N]>>[C:1](=O)N",
    "nitrile → NH2 (reduction)":       "[C:1]#[N]>>[C:1]CN",
    "nitrile → NHCH3":       "[C:1]#[N]>>[C:1]CN(C)",
    "nitrile → N(CH3)2":     "[C:1]#[N]>>[C:1]CN(C)C",
    "nitrile → NHC6H5":      "[C:1]#[N]>>[C:1]Nc1ccccc1",
    "nitrile → C6H4CN":      "[C:1]#[N]>>[C:1]c1ccc(C#N)cc1",
    "nitrile → CHO":         "[C:1]#[N]>>[C:1]=O",

    # 3) Target: Gamma-lactam ring (2-pyrrolidinone)
    "2-pyrrolidinone → 2-piperazinone (template)": "O=C1NCCC1>>O=C1NCCNCC1",

    # 4) Target: -CONH (amide: C(=O)N-)
    "amide → NHCONH2 (urea-like)":      "[C:1](=O)[N:2]>>[C:1](=O)NN",
    "amide → O–C(=O) (carbamate)":      "[C:1](=O)[N:2]>>[C:1](=O)OC",
    "amide → SO2NH2 (sulfonamide-like)": "[C:1](=O)[N:2]>>[C:1]NS(=O)(=O)N",

    "amide → NH2":           "[C:1](=O)[N:2]>>[C:1]N",
    "amide → NHCH3":         "[C:1](=O)[N:2]>>[C:1](=O)N(C)",
    "amide → N(CH3)2":       "[C:1](=O)[N:2]>>[C:1](=O)N(C)C",

    "amide → CHO":           "[C:1](=O)[N:2]>>[C:1](=O)",
    "amide → COOH":          "[C:1](=O)[N:2]>>[C:1](=O)O",
    "amide → CONH2":         "[C:1](=O)[N:2]>>[C:1](=O)N",

    "amide → NHC6H5":        "[C:1](=O)[N:2]>>[C:1](=O)Nc1ccccc1",
    "amide → C6H4CN":        "[C:1](=O)[N:2]>>[C:1](=O)c1ccc(C#N)cc1",

    # 5) Target: -CF3
    "CF3 → SO2NH2":          "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]S(=O)(=O)N",
    "CF3 → CN":              "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]C#N",
    "CF3 → COOH":            "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]C(=O)O",
    "CF3 → CH3":             "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]C",
    "CF3 → OCH3":            "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]OC",
    "CF3 → NH2":             "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]CN",
    "CF3 → CHO":             "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]C=O",
    "CF3 → CONH2":           "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]C(=O)N",
    "CF3 → NHCH3":           "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]CN(C)",
    "CF3 → N(CH3)2":         "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]CN(C)C",
    "CF3 → NHC6H5":          "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]Nc1ccccc1",
    "CF3 → C6H4CN":          "[C:1][C;X3](C(F)(F)F)>>[C:1][C;X3]c1ccc(C#N)cc1",

    # 6) Target: -F
    "F → NH2":               "[C:1]F>>[C:1]N",
    "F → CHO":               "[C:1]F>>[C:1]=O",
    "F → COOH":              "[C:1]F>>[C:1]C(=O)O",
    "F → CONH2":             "[C:1]F>>[C:1]C(=O)N",
    "F → NHCH3":             "[C:1]F>>[C:1]N(C)",
    "F → N(CH3)2":            "[C:1]F>>[C:1]N(C)C",
    "F → C6H4NH2":            "[C:1]F>>[C:1]c1ccc(N)cc1",
    "F → C6H4CN":             "[C:1]F>>[C:1]c1ccc(C#N)cc1",
    "F → SCN":                "[C:1]F>>[C:1]SC#N",
    "F → SC6H5":              "[C:1]F>>[C:1]Sc1ccccc1",
}


# === FUNCTIONS ===
def apply_reactions(mol, round_label, base_name):
    modified_mols = []
    for label, smarts in reaction_smarts_dict.items():
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
            for i, product in enumerate(products):
                mod_mol = product[0]
                Chem.SanitizeMol(mod_mol)
                mol_H = Chem.AddHs(mod_mol)
                AllChem.EmbedMolecule(mol_H, AllChem.ETKDG())
                AllChem.UFFOptimizeMolecule(mol_H)

                name = f"{base_name}__{label.replace(' → ', '_to_').replace(' ', '_')}_{i+1}"
                file_path = os.path.join(round_label, f"{name}.sdf")
                writer = SDWriter(file_path)
                writer.write(mol_H)
                writer.close()

                log_entries.append({
                    "Round": round_label,
                    "Original": base_name,
                    "Transformation": label,
                    "Product_Name": name,
                    "File": file_path
                })
                modified_mols.append((mol_H, name))
        except Exception as e:
            print(f"Failed for {label}: {e}")
    return modified_mols

# === ROUND 1 ===
suppl = Chem.SDMolSupplier(input_sdf)
round1_mols = []
for idx, mol in enumerate(suppl):
    if mol is None:
        continue
    mol_name = f"Mol{idx+1}"
    round1_mols.extend(apply_reactions(mol, round1_dir, mol_name))

# === ROUND 2 ===
for mol, parent_name in round1_mols:
    apply_reactions(mol, round2_dir, parent_name)

# === LOGGING ===
pd.DataFrame(log_entries).to_csv("/mnt/data/transformation_log.csv", index=False)
print("Transformations complete. Logs saved.")
