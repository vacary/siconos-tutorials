

def analytic():
    v1_m = 2

    m_1=1.
    m_2=1e-0

    e=0

    p = - m_1* m_2*(1+e)/(m_1+ m_2)* v1_m
    v2_p = - p /m_2
    v1_p = p /m_1 + v1_m

    print("p = ", p)
    print("v1_p =", v1_p)
    print("v2_p =", v2_p)

    assert(abs(m_1* (v1_p-v1_m) - p ) < 1e-14)
    assert(abs(m_2* (v2_p) + p ) < 1e-14)
    assert(abs((v1_p - v2_p) + e* (v1_m) ) < 1e-14)
    print("\n")


    def percussion(m_1, m_2, e, v1_m):
        return - m_1* m_2*(1+e)/(m_1+ m_2)* v1_m
    for i in range(10):
        m =1./(10**i)
        print("m, p", m, percussion(m_1,m,e,v1_m))

#analytic()


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
