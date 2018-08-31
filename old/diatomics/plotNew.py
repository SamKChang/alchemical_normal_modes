#!/usr/bin/python

import qctoolkit as qtk
import glob, os, re
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
import copy

qtk.setting.quiet = True

out_files = sorted(glob.glob('toLocal/*/*-d*.out'))
out_files_free = sorted(glob.glob('toLocal/*/*-free*.out'))
out_t_files = sorted(glob.glob('toLocal_s3/*/*-d*.out'))

def freeList(sys):
    if sys != 'N2':
        element_strs = re.findall('[A-Z][^A-Z]*', sys)
        elements = [qtk.element[e] for e in element_strs]
    else:
        elements = [qtk.element['N'], qtk.element['N']]
    outs = []
    for e in elements:
        for out in frees:
            if out.element == e:
                outs.append(out)
        #outs.append([out for out in frees if out.element.symbol == e.symbol][0])
    return outs

def getOut(sys, mols, x):
    out_list = [out for out in mols if sys in out.name]
    name_list = [out.name for out in out_list]
    flag = 'd%5.3f.out' % x
    match_list = filter(lambda x: flag in x, name_list)
    if len(match_list) > 0:
        name = match_list[-1]
        ind = name_list.index(name)
        return out_list[ind]
    else:
        return None

def Eele(sys, x_list, mols):
    elements = freeList(sys)
    Si_e = freeList('Si')[0]
    Si_E = Si_e.Et
    freeE = sum([(e.Et) for e in elements])
    if sys != 'Si':
        EList = []
        for x in x_list:
            out = getOut(sys, mols, x)
            if out is not None:
                EList.append(out.Et - out.nuclear_repulsion)
            else:
                EList.append(np.nan)
        # set first energy as silicon
        EList[0] = Si_E
    else:
        EList = [Si_E for i in x_list]

    return EList

def Ediff(sys, x_list, mols1, mols2):
    EList = []
    for x in x_list:
        out1 = getOut(sys, mols1, x)
        out2 = getOut(sys, mols2, x)
        if out1 is not None and out2 is not None:
            diff = out1.Et - out2.Et
            if abs(diff) < 0.5:
                EList.append(diff)
            else:
                EList.append(np.nan)
        else:
            EList.append(np.nan)
    return EList

def plot3DLines(E_function, ax, x, *molList):
    y = range(len(sys_plot))
    X, Y = np.meshgrid(x, y)
    cross_section = []
    cross_ind = []

    new_Z = True
    for i in range(len(sys_plot)):
        sys = sys_plot[i]
        y = [i for _ in d_list]
        z = E_function(sys, x, *molList) # <---- construct electronic energy list
        if new_Z:
            new_Z = False
            Z = np.atleast_2d([z])
        else:
            Z = np.vstack([Z, z])
        cross = []
        for c in cross_coord:
            x_np = np.array(x)
            ind = np.argmin(abs(x_np - c))
            if ind not in cross_ind:
                cross_ind.append(ind)
            cross.append(z[ind])
        cross_section.append(cross)

        ax.plot(x, y, z, color='b')

    cross_section = np.array(cross_section).T
    for i in range(len(cross_section)):
        z = cross_section[i]
        x = [d_list[cross_ind[i]] for _ in z]
        y = range(len(z))
        ax.plot(x,y,z, color='r', ls=':')

    ax.contour(X, Y, Z, offset = E_Si, colors='k')
    return ax

sys_list = []
d_list = []
e_list = []

def loadMol(out_files):
    mols = []
    for out in out_files:
        path, name = os.path.split(out)
        sys, root = name.split('-d')
        d = float(root.replace('.out', ''))
        if sys not in sys_list:
            sys_list.append(sys)
        if d not in d_list:
            d_list.append(d)
        mols.append(qtk.QMOut(out, program='nwchem'))
    return mols

mols = loadMol(out_files)
mols_t = loadMol(out_t_files)

frees = []
for out in out_files_free:
    name = out.split('free-')[-1].replace('.out', '')
    e_list.append(qtk.element[name])
    qout = qtk.QMOut(out, program='nwchem')
    qout.element = e_list[-1]
    frees.append(qout)
    
E_Si = freeList('Si')[0].Et

sys_plot = ['N2', 'CO', 'BF', 'NeBe', 'NaLi', 'HeMg', 'AlH', 'Si']
cross_coord = [0.5, 1.5, 2.5]

# electronic energy contour plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax = plot3DLines(Eele, ax, d_list, mols)      # <---- fill ax plot with data
#plotEele(mols_t, ax)   # <---- fill ax plot with data
y_labels = copy.deepcopy(sys_plot)
y_labels[0] = 'N$_2$'
ax.xaxis.set_rotate_label(False) 
ax.zaxis.set_rotate_label(False) 
ax.set_xlim([0,2.5])
plt.xlim([0,2.5])
ax.set_yticklabels(y_labels,rotation=-20)
ax.set_xlabel('$d$ $[\AA]$', rotation = 0, fontsize=15, labelpad=15)
ax.set_zlabel('$E_{ele}$ $[Ha]$', rotation = -90)
plt.autoscale()
ax.view_init(20, 150)
plt.savefig("Eele.png",bbox_inches='tight')

fig = plt.figure('triplet-singlet')
ax = fig.gca(projection='3d')
ax = plot3DLines(Ediff, ax, d_list, mols, mols_t)
plt.savefig("Ediff.png",bbox_inches='tight')
# limit plot range does not work under ipython

plt.show()
