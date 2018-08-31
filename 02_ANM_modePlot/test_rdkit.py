from rdkit import Chem
from rdkit.Chem import Draw

size = (120, 120)  # Smaller figures than the default.

m = Chem.MolFromSmiles('c1ccccc1')
fig = Draw.MolToMPL(m, size=size)

print fig
