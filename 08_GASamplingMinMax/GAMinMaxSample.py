
# coding: utf-8

# In[1]:

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
ccs_BN = [qtk.CCS(mol_base, yml) for yml in sorted(glob('ccs/coronene_BN*.yml'))]
qtk.setting.quiet=False


# In[23]:

def energy_estimate(dZ):
    # normal mode decomposition
    dZes = dZ.dot(ev_data)
    # calculate first/second order estimate
    d1E = dZ.dot(d1E_Z)
    d2E = (dZes ** 2).dot(ew_data) * 0.5
    return d1E + d2E

def penalty_function(ccs_coord, ccs, name, qmsetting_dict={}):

    time_stamp = datetime.now()

    mol_mut = ccs.generate(**ccs_coord)
    mol_mut.name = '%s_%s' % (name, time_stamp.strftime('%m%d%H%M%S%f')[:11])
    mol_mut.name = mol_mut.name + '_' + str(os.getpid())[-3:]
    
    # construct vector of dZ, lenght of 24 for coronene
    mol_mut_dZ = mol_mut.Z[:N_heavy] - mol_base.Z[:N_heavy]
    D2E = energy_estimate(mol_mut_dZ)
    
    return D2E, mol_mut.name


## In[27]:
#
#ccs = ccs_BN[3]
#coord_test = ccs.random()[1]
#print penalty_function(coord_test, ccs, 'test')
#mol_test = ccs.generate(**coord_test)
#print mol_test.Z


## In[28]:
#
#ccs = ccs_BN[0]
#crd1 = ccs.random()[1]
#crd2 = ccs.random()[1]
#print crd1['mutation'], crd2['mutation']
#Z = np.array(ccs.mate(crd1, crd2, 0.05)['mutation'][0])
#Z[Z != 6]


## In[29]:
#
#test_dz = np.loadtxt('../02_ANM_modePlot/data_space/coronene-dz-02-data.txt')
#test_dE = np.loadtxt('../02_ANM_modePlot/data_space/coronene-dE-02-data.txt')
#energy_estimate(test_dz) - test_dE[:,0]


# In[30]:

def get_optimizer(ccs, mode, name=None, max_step=10000):
    if name is None:
        name = os.path.splitext(ccs.__repr__().replace(' ','_'))[0]
    log_file = '%s.db' % name
    penalty_input = [ccs, name, {}]
    optimizer = qtk.optimization.GeneticOptimizer(
        penalty_function,
        penalty_input,
        ccs.random_coord,#genCCSInp,
        ccs.mate,
        20,
        log_file=log_file,
        new_run=True,
        max_step=max_step,
        mode=mode,
    )
    return optimizer


## In[31]:
#
#test_optimzer = get_optimizer(ccs_BN[1], 'minimize', 'test', max_step=5)
#
#
## In[32]:
#
#test_optimzer.run()


# In[33]:

optimizers = []
max_step = 10
for ccs in ccs_BN:
    name = os.path.splitext(ccs.parameter_file)[0].split('_')[-1]
    opt = get_optimizer(ccs, 'minimize', name+'_min', max_step=max_step)
    optimizers.append(opt)
    #opt.run()
    opt = get_optimizer(ccs, 'maximize', name+'_max', max_step=max_step)
    optimizers.append(opt)
    #opt.run()


# In[14]:

try:
    ind = int(sys.argv[1])
    assert 0 <= ind < 24
    try:
        optimizers[ind].run()
    except Exception as err:
        print str(err)
except:
    print "command line arg 0-23"


# In[ ]:

