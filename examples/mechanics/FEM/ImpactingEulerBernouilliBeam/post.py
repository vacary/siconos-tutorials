import numpy as np

res= np.loadtxt('percussion.txt')
#res= np.loadtxt('percussion.txt.lumped')

#print('res', res)


h_list = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
n_list = [10, 50, 100, 500, 1000, 5000, 10000]

# split result

n_h = len(h_list)
n_el = len(n_list)



import matplotlib.pyplot as plt

print(res[0:4, 0],res[0:4, 2])

plt.subplot(211)
cnt = 0
for n in n_list:
    l = 'n = '+ str(n)
    #print(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3])
    plt.plot(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3], marker='+', label= l )
    cnt =cnt+1
plt.xlabel('h')
plt.ylabel('First percussion')
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.grid()
plt.title('percussion with respect to h')

plt.subplot(212)
cnt = 0
for n in n_list:
    l = 'n = '+ str(n)
    plt.plot(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 4], marker='+', label= l )
    cnt =cnt+1
plt.xlabel('h')
plt.ylabel('W(0,0)')
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.grid()
plt.title('W(0,0) with respect to h')

plt.figure()

plt.subplot(211)
cnt = 0
for n in n_list:
    l = 'n = '+ str(n)
    print('refl value', res[-1,3])
    
    plt.plot(res[cnt*n_h:(cnt+1)*n_h-1, 1], abs(res[cnt*n_h:(cnt+1)*n_h-1, 3]-res[-1,3]), marker='+', label= l )
    cnt =cnt+1
plt.xlabel('h')
plt.ylabel('First percussion')
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.grid()
plt.title('convergence percussion with respect to h')



plt.subplot(212)
cnt = 0
for n in n_list:
    l = 'n = '+ str(n)
    print('refl value', res[-1,4])
    plt.plot(res[cnt*n_h:(cnt+1)*n_h-1, 1], abs(res[cnt*n_h:(cnt+1)*n_h-1, 4]-res[-1,4]), marker='+', label= l )
    cnt =cnt+1
plt.xlabel('h')
plt.ylabel('W(0,0)')
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.grid()
plt.title('convergence W(0,0) with respect to h')



plt.show()
