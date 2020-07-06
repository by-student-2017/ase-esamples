from ase.build import bulk
from ase.calculators.espresso import Espresso
atoms = bulk("Si") #バルクのSiの用意
pseudopotentials = {'Si':'Si.pz-vbc.UPF'} #擬ポテンシャルの設定


input_data = {
    'system': {
        'ecutwfc': 64,
        'ecutrho': 576,
        'nbnd' : 12 },
    'disk_io': 'low'}  #Quantum Espressoのパラメータ

# https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html
calc = Espresso(pseudopotentials=pseudopotentials,kpts=(4, 4, 4),input_data=input_data,pseudo_dir = './')
atoms.set_calculator(calc)

atoms.get_potential_energy()
fermi_level = calc.get_fermi_level()
print(fermi_level)



input_data.update({'calculation':'bands',
                    'restart_mode':'restart',
                    'verbosity':'high'})
calc.set(kpts={'path':'LGXWG', 'npoints':100},
    input_data=input_data)
calc.calculate(atoms)



import matplotlib.pyplot as plt

bs = calc.band_structure()
bs.reference = fermi_level

bs.plot(emax=15, filename='Si.png')

