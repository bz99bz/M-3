import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# Read the CSV file
file_path = './MPP/BACE_multi.csv'
data = pd.read_csv(file_path,encoding='latin1')
data['mol'] = data['mol'].astype(str) 

# Function to convert a SMILES string to an image
def smiles_to_image(smiles, img_path):
    mol = Chem.MolFromSmiles(smiles)
    if mol:  # Check if the molecule is valid
        img = Draw.MolToImage(mol)
        img.save(img_path)
    else:
        print(f"Invalid SMILES string: {smiles}")

# Loop through the DataFrame and save each SMILES as an image
for smiles in data['mol']:
    print(smiles)
    img_path = f"molecule_{smiles}.png"
    smiles_to_image(smiles, img_path)
    print(f"Saved {img_path}")

print("All images have been saved.")
