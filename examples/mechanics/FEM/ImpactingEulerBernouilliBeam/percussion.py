

import subprocess
import os

if os.path.exists('percussion.txt'):
    os.rename('percussion.txt', 'percussion.txt.bck')
    
cmd = ['siconos', 'ImpactingEBBeam.cpp']
returned_value = subprocess.call(cmd, shell=False)  # returns the exit code in unix
print('returned value:', returned_value)
if os.path.exists('percussion.txt'):
    os.remove('percussion.txt')


h_list = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
n_el = [10, 50, 100, 500, 1000, 5000, 10000]
for n in n_el:
    for h in h_list:
        cmd = ['ImpactingEBBeam', str(h), str(n)]
        returned_value = subprocess.call(cmd, shell=False)  # returns the exit code in unix
        print('returned value:', returned_value)
