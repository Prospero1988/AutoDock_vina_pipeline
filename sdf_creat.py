from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp

# CID (PubChem IDs) dla retinolu, kapsaicyny i cholesterolu
compounds = {
    "retinol": 445354,
    "capsaicin": 1548943,
    "cholesterol": 5997
}

# Nazwa pliku wyjściowego
output_sdf = "lig_test.sdf"

# Inicjalizacja writer'a do pliku SDF
writer = Chem.SDWriter(output_sdf)

for name, cid in compounds.items():
    try:
        # Pobranie SMILES z PubChem
        compound = pcp.Compound.from_cid(cid)
        smiles = compound.isomeric_smiles
        if not smiles:
            print(f"Nie udało się pobrać SMILES dla: {name}")
            continue

        # Generowanie molekuły z SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Nie udało się wygenerować molekuły dla: {name}")
            continue

        # Ustawienie współrzędnych 3D
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)

        # Dodanie nazwy jako właściwość w pliku SDF
        mol.SetProp("_Name", name)

        # Zapis do pliku SDF
        writer.write(mol)
        print(f"Struktura dla {name} została dodana do {output_sdf}")

    except Exception as e:
        print(f"Błąd dla {name}: {e}")

# Zamknięcie writer'a
writer.close()

print(f"Plik {output_sdf} został wygenerowany.")
