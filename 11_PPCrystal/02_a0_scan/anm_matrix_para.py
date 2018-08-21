
# coding: utf-8

# # Construct of general ANM matrix using PP
# using PP for ANM is quite a lot more complicated and tedious than all-electron systems.
# This note book attempt to create a readable/reusable script to construct PP ANM matrix

# In[1]:

import qctoolkit as qtk
import numpy as np
import os


# In[2]:

# reuse code from PRM paper
cyl_map = [
        ['si', 5.431, [14, 14]],
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

]
print len(cyl_map)


# In[3]:

# only CPMD PP interface is ready at the moment
# only Gamma-point calculations are required
qmsetting = {
    'program': 'cpmd',
    'kmesh': [1,1,1],
    'link_dep': False,
}


# In[4]:

mol_base = qtk.Molecule('xyz/ge.xyz')
# mol_base.extend([2,2,2])
# mol_base.sort(order='xyz')


# # ANM space and finite difference matrix
# ANM space is spanned by the PP parameters. For FCC primitive cell of two atoms, there are 30 parameters (15 per atoms).
# 
# The goal is the construct matrix elements 
# $$
# \begin{array}{rcl}
# \mathbf{H}_{ij}&=&\frac{\partial^2 E}{\partial\sigma_i\partial\sigma_j}\\
# &=&\frac{\partial}{\partial\sigma_i}\Big(\frac{\partial E}{\partial\sigma_j}\Big)\\
# &\approx&
# \frac{\partial}{\partial\sigma_i}\Big(\frac{E(\sigma_i, \sigma_j+\Delta\sigma_j) - E(\sigma_i, \sigma_j)}{\Delta\sigma_j}\Big)\\
# &\approx&
# \frac{1}{\Delta\sigma_i}\Bigg(\frac{E(\sigma_i+\Delta\sigma_i, \sigma_j+\Delta\sigma_j) - E(\sigma_i+\Delta\sigma_i, \sigma_j)}{\Delta\sigma_j}
# - \frac{E(\sigma_i, \sigma_j+\Delta\sigma_j) - E(\sigma_i, \sigma_j)}{\Delta\sigma_j}\Bigg).
# \end{array}
# $$
# 
# That is, there are four finite difference calculations required for each of the matrix elements:
# $E(\sigma_i+\Delta\sigma_i, \sigma_j+\Delta\sigma_j)$, 
# $E(\sigma_i+\Delta\sigma_i, \sigma_j)$, 
# $E(\sigma_i, \sigma_j+\Delta\sigma_j)$, 
# $E(\sigma_i, \sigma_j)$, 
# where only the first term is unique for each element.
# 
# Note that the finite difference formula is different for diagonal terms $\mathbf{H}_{ii} = \frac{E(\sigma_i+\Delta\sigma_i) - 2E(\sigma_i) + E(\sigma_i-\Delta\sigma_i)}{\Delta\sigma_i^2}$.
# 
# The required finite difference calculations are
# * $E(\sigma_i+\Delta\sigma_i, \sigma_j+\Delta\sigma_j)$: $N(N-1)/2$ calculations for $i\neq j$
# * $E(\sigma_i+\Delta\sigma_i)$: $N$ calculations
# * $E(\sigma_i-\Delta\sigma_i)$: $N$ calculations
# * $E(\sigma_i, \sigma_j)$: 1 calculation
# 
# which added up to $\frac{N^2}{2}+\frac{3}{2}N+1$ calculations, where $N$ is the number of parameters in the system.

# # Function interface
# The following functions are required for reusable/extensive code
# 
# ## _crystal to parameter space map_: cyl2par
# input: crystal object
# 
# output: parameter space coordinate
# 
# ## _parameter spsce to crystal map_: par2cyl
# input: parameter space coordinate
# 
# output: crystal object
# 
# ## _finite difference crystal construction_: FDcyl
# input: finite difference indices, step size
# 
# output: QMInp calculation object
# 
# __NOTE__: for the construction of FD matrix, the index information will be written on the file name for lazy processing
# 
# ## _finite difference matrix construction_: H_FD
# input: finite difference calculation folder path __with fomated index file name__
# 
# output: finite difference matrix as 2D numpy array
# 
# ### TODO:
# implement CPMD interface to take PP object as input and write PP file on the fly __DONE__

# In[5]:

def cyl2par(mol, size=[1, 3]):
    par = []
    for a in mol.type_list:
        crd = qtk.PP(a, size=size)
        vec = crd.vectorize()
        par.append([vec[0][1:], a, vec[0][0]])
    return par

def par2cyl(par, mol_base=mol_base, size=[1, 3]):
    mol = mol_base.copy()
    mol.sort(order='xyz')
    for i, crd in enumerate(par):
        pp = qtk.PP(crd[1], size=size)
        vec = [crd[2]]
        vec.extend(crd[0])
        pp.unvectorize(vec, size[0], range(size[1], 0, -1))
        pp.name = crd[1] + '%02d' % i
        mol.setAtoms(i, string=pp)
    return mol

mol = par2cyl(cyl2par(mol_base))
inp = qtk.QMInp(mol, **qmsetting)
inp.write()


# In[6]:

def FDcyl(mol, i, j, step, size=[1,3], dim=15):
    """
    input:  molecule object: where the FD perturbation will be introduced
            i, j: integer indices for the perturbation matrix row and column number
            step: absolute step size for finite difference perturbation
            size: PP configuration space
    output: molecule object with appropreate PP and file names attached
    """
    def get_pp(I, ip, font, mode='+'):
        if type(mol.string[I]) is type(qtk.PP()):
            pp = mol.string[I]
        else:
            pp = qtk.PP(mol.type_list[I], size=size)
            pp.name = mol.type_list[I]
        vec = pp.vectorize()[0]
        new_pp = pp.copy()
        if mode == '+':
            vec[1+ip] += step
        elif mode == '-':
            vec[1+ip] -= step
        new_pp.unvectorize(vec, size[0], range(size[1], 0, -1))
        new_pp.name += '_%s%02dp%02d' % (font, I, ip)
        return new_pp
    
    # first order
    if i == j:
        I, ip = divmod(i, dim)
        modes = [['p', '+'], ['m', '-']]
        mols = []
        for m in modes:
            mol_new = mol.copy()
            pp = get_pp(I, ip, 'I', mode=m[1])
            mol_new.setAtoms(I, string=pp)
            mol_new.name += '_i%02d%s-I%02dp%02d' % (i, m[0], I, ip)
            mols.append(mol_new)
        return mols
    # second order
    else:
        mol = mol.copy()
        inds = [i, j]
        ind_font = ['I', 'J']
        for t, s in enumerate([i, j]):
            I, ip = divmod(s, dim)
            inds.append(I)
            inds.append(ip)
            pp = get_pp(I, ip, ind_font[t])
            mol.setAtoms(I, string=pp)
        mol.name += "_i%02dj%02d-I%02dp%02d-J%02dp%02d" % tuple(inds)
    return [mol]

# mol_test = FDcyl(mol_base, 5, 28, 100)[0]
# inp = qtk.QMInp(mol_test, **qmsetting)
# inp.write('PPTest', overwrite=True)


# ## Construct FD crystal for calculations

# In[ ]:

mols_fd = []
par = cyl2par(mol_base)

# get dimension of PP space
# single atom PP space * number of atoms
N = len(par[0][0]) * len(par)

for info_cyl in cyl_map:
    mols_c = []
    mols_info = [info_cyl[0] + '_H_scan', mols_c]
    if not os.path.exists(mols_info[0]):
        mols_fd.append(mols_info)
        print "processing %s" % info_cyl[0]
        mol_c = mol_base.copy()
        for i in range(mol.N):
            mol_c.setAtoms(i, Z=info_cyl[2][i])
        mol_c.name = info_cyl[0]
        for info_a in cyl_map:
            mol_a = mol_c.copy()
            mol_a.name += '_a-' + info_a[0]
            mol_a.setCelldm([info_a[1], info_a[1], info_a[1], 0, 0, 0])
            mols_c.append(mol_a)
            # get off-diagonal terms
            for i in range(N):
                for j in range(i, N):
                    mols = FDcyl(mol_a, i, j, step=0.01)
                    mols_c.extend(mols)
        inps_c = [qtk.QMInp(m, **qmsetting) for m in mols_c]
        qtk.qmRunAll(inps_c, mols_info[0])
    else:
        print "%s exists, skip" % info_cyl[0]
        



