#!/usr/bin/env python

from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS
import ase.io 

pseudopotentials = {'Na': 'Na.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                            'Cl': 'Cl.pbesol-n-kjpaw_psl.1.0.0.UPF'}  
rocksalt = bulk('NaCl', crystalstructure='rocksalt', a=6.0)
calc = Espresso(pseudopotentials=pseudopotentials,pseudo_dir = './',
                        tstress=True, tprnfor=True, kpts=(3, 3, 3))

rocksalt.set_calculator(calc)

ucf = UnitCellFilter(rocksalt)
opt = LBFGS(ucf)
opt.run(fmax=0.005)

# cubic lattic constant
print((8*rocksalt.get_volume()/len(rocksalt))**(1.0/3.0))
