
# coding: utf-8

# In[2]:

import qctoolkit as qtk
import numpy as np
import glob
import os, re
import copy


# In[3]:

cyl_map = [
        ['si', 5.431, [14, 14]],
#        ['c', 3.566, [6, 6]],
#        ['sic', 4.3596, [14, 6]],
#        ['cge', 5.574, [6, 32]],
        ['ge', 5.658, [32, 32]],
        ['sige', 5.432, [14, 32]],
        ['sn', 6.4892, [50, 50]],
        ['snsi', 5.961, [14, 50]],
        ['gesn', 6.0758, [32, 50]],
        ['gaas', 5.6535, [31, 33]],
        ['alp', 5.4635, [13, 15]],
        ['alas', 5.660, [13, 33]],
        ['alsb', 6.1355, [13, 51]],
        ['gap', 5.451, [31, 15]],
        ['gasb', 6.09, [31, 51]],
        ['inp', 5.86, [49, 15]],
        ['inas', 6.05, [49, 33]],
        ['insb', 6.47, [49, 51]],
# II-VI has different number of valence electrons
#         ['cdse', 6.08, [48, 34]],
#         ['zns', 5.41, [30, 16]],
#         ['znse', 5.67, [30, 34]],
#         ['znte', 6.1, [30, 52]],
#         ['cds', 5.82, [48,16]],
#         ['cdte', 6.48, [48, 52]],
#         ['hgs', 5.8517, [80, 16]],
#         ['hgse', 6.085, [80, 34]],
#         ['hgte', 6.453, [80, 52]],
]
print len(cyl_map)


# # plain GGA, where insuite restart is possible

# In[4]:

# # single setting dict for insuit restart calculations
# qmsetting = {
#     'program': 'abinit',
#     'kmesh': [6,6,6],
#     'ks_states': 10,
#     'band_scan': [
#      [17, 17, 8, 10, 15],
#      [[0.25, 0.75, 0.5], # W
#       [0.0, 0.0, 0.0],   # Gamma
#       [0.0, 0.5, 0.5],   # X
#       [0.25, 0.75, 0.5], # L
#       [0.5, 0.5, 0.5],   # W
#       [0.0, 0.0, 0.0],   # Gamma
#      ]],
#     'save_restart': True,
#     'save_density': True,
#     'link_dep': True,
#     'threads': 1,
#     'dos_mesh': [12, 12, 12],
#     'smearing': 'gaussian',
# }

qmsetting = {
    'program': 'cpmd',
    'kmesh': [1,1,1],
    'link_dep': True,
}


# In[5]:

# construct crystals based on cyl_map and si, gaas, cdse crystals

mol_base = qtk.Molecule('xyz/gaas.xyz')


# In[17]:

mols = []
count = 0
for info_mol in cyl_map:
    print info_mol[0]
    count += 1
    for info_a in cyl_map:
        new = mol_base.copy()
        new.name = info_mol[0] + '_a-' + info_a[0]
        new.setCelldm([info_a[1], info_a[1], info_a[1], 0, 0, 0])
        for i in range(new.N):
            if new.Z[i] != info_mol[2][i]:
                new.setAtoms(i, Z=info_mol[2][i])
        mols.append(new)
print count
#mols_all = list(qtk.flatten(mols_grp))
inps_ref = [qtk.QMInp(mol, **qmsetting) for mol in mols]
print len(inps_ref)
print inps_ref[20]


# In[14]:

# # uniform PP space
# Zs = [13, 14, 15, 31, 32, 33, 49, 50, 51]
# for Z in Zs:
#     elem = qtk.Z2n(Z)
#     pp = qtk.PP(elem).resize(1, [3,2,1])
#     print elem, pp,

