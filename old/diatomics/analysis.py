#!/usr/bin/python

import qctoolkit as qtk
import glob, os

outs = sorted(glob.glob('*/*.out'))

sys_list = []
d_list = []
mols = []
for out in outs:
  path, name = os.path.split(out)
  sys, root = name.split('-d')
  d = float(root.replace('.out', ''))
  if sys not in sys_list:
    sys_list.append(sys)
  if d not in d_list:
    d_list.append(d)
  mols.append(qtk.QMOut(out, program='nwchem'))

def EtSys(sys):
  Et = [mol.Et for mol in mols if sys in mol.name]
  return Et

print EtSys(sys_list[0])
