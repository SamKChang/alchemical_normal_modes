#!/usr/bin/python

import qctoolkit as qtk
import numpy as np

qmsetting = {
    'program': 'nwchem',
    'basis_set': 'aug-cc-pVTZ',
    'theory': 'pbe0',
    'threads': 1,
}

def freeMol(A_list, name_root, i):
  R = [[0,0,0]]
  Z = A_list[i].number
  if Z > 0:
    mol = qtk.Molecule()
    mol.name = name_root + '-free-' + A_list[i].symbol
    mol.addAtoms(Z, R)
    if mol.getValenceElectrons() % 2 != 0: 
      mol.setChargeMultiplicity(0, 2)
    mol.write()
    return mol

Z1 = range(7, 15)
Z2 = range(8)
Z2.reverse()

A1 = [qtk.element[z] for z in Z1]
A2 = [qtk.element[z] for z in Z2]

d_list = []
n_list = []
for i in range(len(A1)):
    a = A1[i]
    b = A2[i]
    if b.number > 0:
      if a.eleaffin < b.eleaffin:
          name = a.symbol + b.symbol
      elif b.eleaffin < a.eleaffin:
          name = b.symbol + a.symbol
      else:
          if b.eleneg < a.eleneg:
              name = b.symbol + a.symbol
          elif a.eleneg < b.eleneg:
              name = a.symbol + b.symbol
          else:
              name = a.symbol + str(2)
    else:
      name = a.symbol
    n_list.append(name)
    d_list.append(a.covrad + b.covrad)

inps = []
for i in range(len(n_list)):
    name_root = n_list[i]
    if A2[i].number > 0:
        for d in np.linspace(0.6, 2, 20):
            l = d_list[0]*d
        
            R = [[0, 0, 0], [l, 0, 0]]
            Z = [a.number for a in [A1[i], A2[i]]]
            mol = qtk.Molecule()
            mol.name = name_root + '-d%5.3f' % l
            mol.addAtoms(Z, R)
            inps.append(qtk.QMInp(mol, **qmsetting))
        mol_free1 = freeMol(A1, name_root, i)
        mol_free2 = freeMol(A2, name_root, i)
        mol_free1.write()
        inps.append(qtk.QMInp(mol_free1, **qmsetting))
        inps.append(qtk.QMInp(mol_free2, **qmsetting))
    else:
        R = [[0,0,0]]
        Z = A1[i].number
        mol = qtk.Molecule()
        mol.name = name_root + '-d%5.3f' % 0
        mol.addAtoms(Z, R)
        mol_free1 = freeMol(A1, name_root, i)
        inps.append(qtk.QMInp(mol_free1, **qmsetting))

        
jobs = []
for inp in inps:
    #print inp.setting
    inp.write()
    jobs.append([inp, inp.molecule.name])

qtk.parallelize(qtk.qmRunJob, jobs)
