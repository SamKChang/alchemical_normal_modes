
# coding: utf-8

# In[11]:

import qctoolkit as qtk
import numpy as np
from glob import glob
from datetime import datetime
import os
import sys


# In[2]:

mol_name = 'coronene'
mol_base = qtk.Molecule(mol_name + '.xyz')

mol_data = np.load('data_' + mol_name + '.npz')
d1E_data = np.load('data_' + mol_name + '_1st.npz')
ev_data = mol_data['v']
ew_data = mol_data['w']
d1E_Z = d1E_data['arr_0']
N_heavy = int(np.ones(mol_base.N)[mol_base.Z > 1].sum())

qtk.setting.quiet=True
ccs_BN = [qtk.CCS(mol_base, yml) for yml in sorted(glob('coronene_BN*.yml'))]
qtk.setting.quiet=False


# In[3]:

def genCCSInp():
    """rapper for GA optimizer"""
    _coord = ccs.random()[1]
    return _coord

#penalty_input = [ccs, {}]

def energy_estimate(dZ):
    # normal mode decomposition
    dZes = dZ.dot(ev_data)
    # calculate first/second order estimate
    d1E = dZ.dot(d1E_Z)
    d2E = (dZes ** 2).dot(ew_data) * 0.5
    return d1E + d2E

def penalty_function(ccs_coord, ccs, target, name, qmsetting_dict={}):

    time_stamp = datetime.now()

    mol_mut = ccs.generate(**ccs_coord)
    mol_mut.name = '%s_%s' % (name, time_stamp.strftime('%m%d%H%M%S%f')[:11])
    mol_mut.name = mol_mut.name + '_' + str(os.getpid())[-3:]
    
    # construct vector of dZ, lenght of 24 for coronene
    mol_mut_dZ = mol_mut.Z[:N_heavy] - mol_base.Z[:N_heavy]
    D2E = energy_estimate(mol_mut_dZ)
    
    return np.sqrt((D2E - target) ** 2), mol_mut.name


## In[4]:
#
#ccs = ccs_BN[3]
#coord_test = ccs.random()[1]
#print penalty_function(coord_test, ccs, -1.1, 'test')
#mol_test = ccs.generate(**coord_test)
#print mol_test.Z


# In[5]:

#test_dz = np.loadtxt('../02_ANM_modePlot/data_space/coronene-dz-02-data.txt')
#test_dE = np.loadtxt('../02_ANM_modePlot/data_space/coronene-dE-02-data.txt')
#energy_estimate(test_dz) - test_dE[:,0]


# In[6]:

def get_optimizer(ccs, target, name=None, max_step=10000):
    if name is None:
        name = os.path.splitext(ccs.__repr__().replace(' ','_'))[0]
    log_file = '%s.db' % name
    penalty_input = [ccs, target, name, {}]
    optimizer = qtk.optimization.GeneticOptimizer(
        penalty_function,
        penalty_input,
        genCCSInp,
        ccs.mate,
        20,
        log=log_file,
        new_run=True,
        max_step=max_step,
    )
    return optimizer


# In[7]:

#test_optimzer = get_optimizer(ccs_BN[3], -1.1, 'test', max_step=10)
#test_optimzer.run()


# In[8]:

optimizers = []
max_step = 10000
for ccs in ccs_BN:
    name = os.path.splitext(ccs.parameter_file)[0].split('_')[-1]
    opt_min = get_optimizer(ccs, -100, name+'_min', max_step=max_step)
    opt_max = get_optimizer(ccs, -100, name+'_max', max_step=max_step)
    optimizers.append(opt_min)
    optimizers.append(opt_max)

for opt in optimizers:
  opt.run()

## In[14]:
#
#try:
#    ind = int(sys.argv[1])
#    assert 0 <= ind < 24
#    try:
#        opt = optimizers[ind]
#        print "runing %s" % opt.log
#        optimizers[ind].run()
#    except Exception as err:
#        print str(err)
#except:
#    print "command line arg 0-23"
#
#
## In[ ]:



