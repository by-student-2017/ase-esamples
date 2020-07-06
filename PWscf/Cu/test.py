#!/usr/bin/env python

import matplotlib.pyplot as plt
from ase import Atoms
from ase.build import bulk
from ase.calculators.espresso import Espresso

atoms = bulk("Cu")
pseudopotentials = {'Cu':'Cu.pz-d-rrkjus.UPF'}

input_data = {
    'system': {
        'ecutwfc': 30,
        'ecutrho': 240,
        'nbnd' : 35,
        'occupations' : 'smearing',
        'smearing':'gauss',
        'degauss' : 0.01},
    'disk_io': 'low'}  # automatically put into 'control'

# https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html
calc = Espresso(pseudopotentials=pseudopotentials,kpts=(4, 4, 4),input_data=input_data,pseudo_dir = './')
atoms.set_calculator(calc)

atoms.get_potential_energy()
fermi_level = calc.get_fermi_level()
print(fermi_level)

input_data.update({'calculation':'bands',
                        'restart_mode':'restart',
                        'verbosity':'high'})
calc.set(kpts={'path':'GXWLGK', 'npoints':100},
    input_data=input_data)
calc.calculate(atoms)

bs = calc.band_structure()
bs.reference = fermi_level
bs.plot(emax=40,emin=5)

