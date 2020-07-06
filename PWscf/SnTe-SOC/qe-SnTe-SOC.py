#!/usr/bin/env python

from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.constraints import UnitCellFilter
import ase.io 

pseudopotentials = {'Sn': 'Sn.rel-pbesol-dn-rrkjus_psl.0.2.UPF',
                    'Te': 'Te.rel-pbesol-dn-rrkjus_psl.0.2.2.UPF'}

atoms =ase.io.read("SnTe_mp-1883_computed.cif")

input_data = {
    'system': {
        'ecutwfc': 30,
        'ecutrho': 240,
        'nbnd' : 35,
        'occupations' : 'smearing',
        'smearing':'gauss',
        'degauss' : 0.01,
        'noncolin': True,
        'lspinorb': True},
    'CONTROL':{
        'calculation':'scf',
        'prefix':'SnTe' ,
        'outdir':'./',
        'pseudo_dir':'./'},
        'disk_io': 'low'} 

calc = Espresso(pseudopotentials=pseudopotentials,kpts=(4, 4, 4),input_data=input_data,pseudo_dir = './')
atoms.set_calculator(calc)


atoms.get_potential_energy()
fermi_level = calc.get_fermi_level()
print(fermi_level)


input_data.update({'calculation':'bands',
                    'restart_mode':'restart',
                    'verbosity':'high'})
calc.set(kpts={'path':'GLUWLU', 'npoints':100},
            input_data=input_data)
calc.calculate(atoms)



import matplotlib.pyplot as plt

bs = calc.band_structure()
bs.reference = fermi_level

bs.plot(emax=13,emin=4,filename='SnTe_SO.png')
bs.plot(emax=fermi_level+2,emin=fermi_level-2,filename='SnTe_mag_SO.png')
