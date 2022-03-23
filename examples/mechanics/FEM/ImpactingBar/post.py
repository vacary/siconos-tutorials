import numpy as np

res= np.loadtxt('percussion.txt')
#res= np.loadtxt('percussion.txt.lumped')
#res= np.loadtxt('percussion_50000.txt')
#res= np.loadtxt('percussion_100000.txt')
#res= np.loadtxt('percussion_new.txt')
#print('res', res)

h_list = list(set(res[:,1].tolist()))
h_list.sort(reverse=True)
print('h_list', h_list)

n_list = list(set((res[:,0]).astype(int).tolist()))
n_list.sort()
print('n_list', n_list)


n_h = len(h_list)
n_el = len(n_list)

# recompute reaction force

res_extended =np.zeros((res.shape[0], res.shape[1]+1))

import math

res_extended[:,:-1] = res[:,:] 
for r in range(res.shape[0]):
    print('row', r, res[r,:])

    res_extended[r,5] = res_extended[r,3] / math.sqrt(res_extended[r,1])
    res_extended[r,5] = res_extended[r,3] / (res_extended[r,1])
    print('row_extended', r, res_extended[r,:])

    
res=res_extended

#input()
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

############# Convergence in h

ax = plt.subplot(211)
# ax.yaxis.set_major_locator(LogLocator(base=100))
# ax.xaxis.set_major_locator(LogLocator(base=100))
# ax.yaxis.set_minor_locator(LogLocator(base=100))
# ax.xaxis.set_minor_locator(LogLocator(base=100))
cnt = 0
for n in n_list:
    l = 'n = '+ str(n)
    
    print(res[cnt*n_h:(cnt+1)*n_h, 0], res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3])
    plt.plot(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3], marker='+', label= l )
    cnt =cnt+1
plt.xlabel('h')
plt.ylabel('First percussion')
plt.yscale("log")
plt.xscale("log")
plt.axis('equal')

plt.legend()
plt.grid()
plt.title('percussion with respect to h')

plt.subplot(212)
cnt = 0
for n in n_list:
    l = 'n = '+ str(n)
    plt.plot(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 5], marker='+', label= l )
    cnt =cnt+1
plt.xlabel('h')
plt.ylabel('reaction p/h')
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.grid()
plt.title('reaction with respect to h')

# plt.subplot(212)
# cnt = 0
# for n in n_list:
#     l = 'n = '+ str(n)
#     plt.plot(res[cnt*n_h:(cnt+1)*n_h, 4], res[cnt*n_h:(cnt+1)*n_h, 4], marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('h')
# plt.ylabel('W(0,0)')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('W(0,0) with respect to h')


# plt.figure()
# cnt=len(n_list)-1
# offset=6
# print(res[cnt*n_h:(cnt+1)*n_h-offset, 1], res[cnt*n_h:(cnt+1)*n_h-offset, 5])
# plt.plot(res[cnt*n_h:(cnt+1)*n_h-offset, 1], res[cnt*n_h:(cnt+1)*n_h-offset, 5], marker='+', label= l )
# plt.xlabel('h')
# plt.ylabel('reaction p/h')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()





# plt.figure()
# plt.subplot(211)
# cnt = 0
# for n in n_list:
#     l = 'n = '+ str(n)
#     print('ref value for n ', res[(cnt+1)*n_h-1,0], ' and h ', res[(cnt+1)*n_h-1,1], ':', res[(cnt+1)*n_h-1,3])
#     print(res[cnt*n_h:(cnt+1)*n_h-1, 3])
#     plt.plot(res[cnt*n_h:(cnt+1)*n_h-1, 1], abs((res[cnt*n_h:(cnt+1)*n_h-1, 3]-res[(cnt+1)*n_h-1,3]))/abs((res[(cnt+1)*n_h-1,3])), marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('h')
# plt.ylabel('First percussion')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('convergence percussion with respect to h')



# plt.subplot(212)
# cnt = 0
# for n in n_list:
#     l = 'n = '+ str(n)
#     print('ref value for n ', res[(cnt+1)*n_h-1,0], ' and h ', res[(cnt+1)*n_h-1,1], ':', res[(cnt+1)*n_h-1,5])
    
#     print(res[cnt*n_h:(cnt+1)*n_h-1, 1], abs((res[cnt*n_h:(cnt+1)*n_h-1, 5]-res[(cnt+1)*n_h-1,5])/abs(res[(cnt+1)*n_h-1,5])))
    
#     plt.plot(res[cnt*n_h:(cnt+1)*n_h-1, 1], abs((res[cnt*n_h:(cnt+1)*n_h-1, 5]-res[(cnt+1)*n_h-1,5])/abs(res[(cnt+1)*n_h-1,5])), marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('h')
# plt.ylabel('reaction p/h')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('convergence reaction with respect to h')

# plt.subplot(212)
# cnt = 0
# for n in n_list:
#     l = 'n = '+ str(n)
#     print('ref value for n ', res[(cnt+1)*n_h-1,0], ' and h ', res[(cnt+1)*n_h-1,1], ':', res[(cnt+1)*n_h-1,5])
#     print(res[cnt*n_h:(cnt+1)*n_h-1, 5])
#     plt.plot(res[cnt*n_h:(cnt+1)*n_h-1, 1], abs((res[cnt*n_h:(cnt+1)*n_h-1, 5]-res[(cnt+1)*n_h-1,5])/abs(res[(cnt+1)*n_h-1,5])), marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('h')
# plt.ylabel('W(0,0)')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('convergence W(0,0) with respect to h')



# ############# Convergence in n


# print('\n ############# Convergence in n \n ')
# # swap result

# res_ori =res.copy()

# res[:,0] = res_ori[:,1]
# res[:,1] = res_ori[:,0]

# print('res', res)
# res_copy =res.copy()

# cnt =0 
# for h in h_list:
#     print('h', h )
#     print( res[cnt::n_h, :])
#     res_copy[cnt*n_el:(cnt+1)*n_el, :] = res[cnt::n_h, :]
#     cnt=cnt+1
#     #input()

# res=res_copy
# print('res', res)



# plt.figure()

# plt.subplot(211)
# cnt = 0
# for h in h_list:
#     l = 'h = '+ str(h)
#     #print(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3])
#     plt.plot(res[cnt*n_el:(cnt+1)*n_el, 1], res[cnt*n_el:(cnt+1)*n_el, 3], marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('n_el')
# plt.ylabel('First percussion')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('percussion with respect to n')


# plt.subplot(212)
# cnt = 0
# for h in h_list:
#     l = 'h = '+ str(h)
#     #print(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3])
#     plt.plot(res[cnt*n_el:(cnt+1)*n_el, 1], res[cnt*n_el:(cnt+1)*n_el, 3], marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('n_el')
# plt.ylabel('W(0,0)')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('W(0,0) with respect to n')




# plt.figure()

# plt.subplot(211)
# cnt = 0
# for h in h_list:
#     l = 'h = '+ str(h)
#     #print(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3])
#     print('ref value for h ', res[(cnt+1)*n_el-1,0], ' and n ', res[(cnt+1)*n_el-1,1], ':', res[(cnt+1)*n_el-1,3])
#     print(abs(res[cnt*n_el:(cnt+1)*n_el-1, 3]-res[(cnt+1)*n_el-1,3]))
#     print(res[cnt*n_el:(cnt+1)*n_el, :])
#     print(res[cnt*n_el:(cnt+1)*n_el-1, 3]-res[(cnt+1)*n_el-1,3])
#     plt.plot(res[cnt*n_el:(cnt+1)*n_el-1, 1], abs(res[cnt*n_el:(cnt+1)*n_el-1, 3]-res[(cnt+1)*n_el-1,3])/abs(res[(cnt+1)*n_el-1,3]), marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('n_el')
# plt.ylabel('Convergence First percussion')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('convergence percussion with respect to n')

# plt.subplot(212)
# cnt = 0
# for h in h_list:
#     l = 'h = '+ str(h)
#     #print(res[cnt*n_h:(cnt+1)*n_h, 1], res[cnt*n_h:(cnt+1)*n_h, 3])
#     print('ref value for h ', res[(cnt+1)*n_el-1,0], ' and h ', res[(cnt+1)*n_el-1,1], ':', res[(cnt+1)*n_el-1,4])
#     print(abs(res[cnt*n_el:(cnt+1)*n_el-1, 4]-res[(cnt+1)*n_el-1,4]))
#     plt.plot(res[cnt*n_el:(cnt+1)*n_el-1, 1], abs(res[cnt*n_el:(cnt+1)*n_el-1, 4]-res[(cnt+1)*n_el-1,4]), marker='+', label= l )
#     cnt =cnt+1
# plt.xlabel('n_el')
# plt.ylabel('Convergence First percussion')
# plt.yscale("log")
# plt.xscale("log")
# plt.legend()
# plt.grid()
# plt.title('convergence W(0,0) with respect to n')


plt.show()
