import qctoolkit as qtk

mol_base = qtk.Molecule('coronene.xyz')
ccs = qtk.CCS(mol_base, 'coronene_BN01.yml')

mol = ccs.random()[0]
print mol.Z
