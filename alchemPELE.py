# Imports
import pathlib 
import os

from peleffy.topology import Molecule 
from peleffy.topology import Topology,Alchemizer
from peleffy.template import Impact
from peleffy.forcefield import OpenForceField
from peleffy.solvent import OBC2

#path = str(pathlib.Path().absolute())
#path_simulation = path + '/HybridPELE'
#path_simulation_inputs = path_simulation + '/input'     
#path_simulation_rotamers = path_simulation + '/rotamers'
#path_simulation_solvent = path_simulation + '/solvent'
#
#if  os.path.exists(path_simulation) == False:
#    os.mkdir(path_simulation)
#elif  os.path.exists(path_simulation_inputs) == False:
#    os.mkdir(path_simulation_inputs)
#elif  os.path.exists(path_simulation_rotamers) == False:
#    os.mkdir(path_simulation_rotamers)
#elif  os.path.exists(path_simulation_solvent) == False:
#    os.mkdir(path_simulation_solvent)

os.environ['SCHRODINGER'] = '/opt/schrodinger2021-3'

# Defining the forcefield we are going to use
openff = OpenForceField('openff_unconstrained-2.0.0.offxml')

# Loading data, parametrizing and topology
m1 = Molecule('4a.pdb')
m2 = Molecule('13k.pdb')

p1 = openff.parameterize(m1, charge_method='gasteiger')
p2 = openff.parameterize(m2, charge_method='gasteiger')

#t1 = Topology(m1, p1)
#t2 = Topology(m2, p2)

#alchemizer = Alchemizer(t1, t2)
#
## Generating files for each step
#for l in range(0, 11):
#    top = alchemizer.get_alchemical_topology
#    impact = Impact(top)
#    n = str(l).rjust(2, '0')
#    impact.to_file(f'4a_4o_{n}')
#    alchemizer.rotamer_library_to_file(f'HYB{n}.rot.assign')
#
## Solvent model
#obc = OBC2(top)
#obc.to_file('ligandParams.txt')
#
#alchemizer.hybrid_to_pdb('hybrid.pdb')