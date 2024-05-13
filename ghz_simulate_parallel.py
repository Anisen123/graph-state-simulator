import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt

def draw_graph(m):
    m = np.matrix(m)
    G = nx.from_numpy_array(m)
    nx.draw_kamada_kawai(G,with_labels = True)

def get_nbs(m,v):
    nbs = []
    for i in range(len(m)):
        if(m[v][i] == 1):
            nbs.append(i)
    return nbs

def flip(x):
    return (x+1)%2

def z(m,v,z_str):
    z_str[v] = flip(z_str[v])

def x(m,v,z_str):
    for i in get_nbs(m,v):
        z(m,i,z_str)
        
def y(m,v,z_str):
    x(m,v,z_str)
    z(m,v,z_str)

def id(m,v,z_str):
    pass

def add_error(m,bit_str,z_str,p_e):
    # If considering memory decoherence as well, uncomment the following 2 lines
    #for i in range(len(m)):
    #    np.random.choice([x,y,z,id],p=[p_e/3,p_e/3,p_e/3,1-p_e])(m,i,z_str)
    for i in bit_str: #(gate noise)
        np.random.choice([x,y,z,id],p=[p_e/3,p_e/3,p_e/3,1-p_e])(m,i,z_str)

#Test add_error
m1 = [[0,1,0,1],[1,0,1,1],[0,1,0,0],[1,1,0,0]]
z_str1 = [0 for i in m1]
add_error(m1,z_str1,z_str1,0.5)
print(z_str1)
        
def cz(m,a,b,z_str,p_e):
    m[a][b] = flip(m[a][b])
    m[b][a] = m[a][b]
    add_error(m,[a,b],z_str,p_e)

def z_meas(m,v,z_str,p_e):
    nbs = get_nbs(m,v)
    for i in nbs:
        m[i][v] = flip(m[v][i])
        m[v][i] = m[i][v]
    add_error(m,[v]+nbs,z_str,p_e)

def y_meas(m,v,z_str,p_e):
    lc(m,v,z_str,p_e)
    nbs = get_nbs(m,v)
    for i in nbs:
        m[i][v] = flip(m[v][i])
        m[v][i] = m[i][v]

def lc(m,v,z_str,p_e):
    nbs = get_nbs(m,v)
    for i in range(len(nbs)):
        for j in range(i+1,len(nbs)):
            m[nbs[i]][nbs[j]] = flip(m[nbs[i]][nbs[j]])
            m[nbs[j]][nbs[i]] = m[nbs[i]][nbs[j]]
    add_error(m,[v]+nbs,z_str,p_e)

############################################################################################################

def par_ghz(m,bit_str,err_str,p_e):
    n = len(bit_str)
    cz(m,bit_str[0],bit_str[1],err_str,p_e)
    k=1
    while 2*k<n:
        for i in range(1,k+1):
            cz(m,bit_str[i],bit_str[k+i],err_str,p_e)
            lc(m,bit_str[i],err_str,p_e)
            cz(m,bit_str[i],bit_str[k+i],err_str,p_e)
        k*=2

    l = k
    if l%2 == 1:
        cz(m,bit_str[0],bit_str[n-1],err_str,p_e)
        l-=1
    for i in range(1,n-k):
        cz(m,bit_str[i],bit_str[l+i],err_str,p_e)
        lc(m,bit_str[i],err_str,p_e)
        cz(m,bit_str[i],bit_str[l+i],err_str,p_e)

    return m

############################################################################################################

#Helper function for SC protocol
def sc_ghz_par(m,n,k,bit_str,err_str,p_e):
    for i in range(1,math.ceil(n//(k-1))+2):
        for j in range(1,k):
            if(j+(k-1)*i>=n+k):
                break
            cz(m,bit_str[j],bit_str[j+(k-1)*i],err_str,p_e)
            lc(m,bit_str[j],err_str,p_e)
            cz(m,bit_str[j],bit_str[j+(k-1)*i],err_str,p_e)


############################################################################################################

def sc(m,n,k,bit_str,err_str,p_e,last = False, reset = True): #m is a matrix of size (n+k)
    assert(k<=n+1 and k>=2)

    #create inner ghz state
    m = par_ghz(m,[i for i in range(k)],err_str,p_e)

    sc_ghz_par(m,n,k,bit_str,err_str,p_e)

    for i in range(k,n+k):
        lc(m,bit_str[i],err_str,p_e)
    
    sc_ghz_par(m,n,k,bit_str,err_str,p_e)

    for i in range(1,k):
        z_meas(m,i,err_str,p_e)
    
    lc(m,bit_str[0],err_str,p_e)

    if(last):
        for i in range(k,n+k):
            z_meas(m,i,err_str,p_e)

        err_str = [err_str[0]]+err_str[k:n+k]
        return int(err_str == [0]*(n+1))

    if(reset):
        #Edge reset operation
        for i in range(1,n):
            cz(m,bit_str[0],bit_str[i],err_str,p_e)

        for i in range(1,n):
            lc(m,bit_str[i],err_str,p_e)

        for i in range(1,n):
            cz(m,bit_str[0],bit_str[i],err_str,p_e)


############################################################################################################
            
#Parallel Subgraph Complementation Algorithm to create GHZ state
def sc_build_ghz_parallel(n,k,p):    
    m1 = [[0 for i in range(2*n+k)] for j in range(2*n+k)]
    err_str = [0 for i in m1]

    for i in range(n):
        cz(m1,i+k,n+k+i,err_str,p) #initialise graph state with bell pairs

    return sc(m1,n,k,list(range(n+k)),err_str,p,True)

############################################################################################################

#Parallel Factory node + Teleportation Protocol to build GHZ state
def cn_build_ghz_parallel(n,p_e):
    m = [[0 for i in range(3*n)] for j in range(3*n)]
    err_str = [0 for i in m]
    
    for i in range(n):
        cz(m,i+n,i+2*n,err_str,p_e) #initialise graph state with bell pairs

    cols = []
    for z in range(1,n+1):
        for i in range(1,n+1):
            j = (z - i) % n
            if j == 0 or j == i:
                continue
            if((i,j) in cols or (j,i) in cols):
                continue
            cz(m,i-1,j-1,err_str,p_e)
            cols.append((i,j))

    for i in range(n):
        cz(m,i,i+n,err_str,p_e)

    for i in range(n):
        y_meas(m,i,err_str,p_e)
        y_meas(m,i+n,err_str,p_e)

    return 1 if err_str[2*n:] == [0]*n else 0

############################################################################################################

# Simulation of parallelized GHZ state generation using SC and CN protocols
# The following code generates graphs of error probability vs fidelity for SC and CN protocols, but this can be easily modified to generate graphs of number of qubits vs fidelity with small changes
sc_arr = []
cn_arr = []
p_arr = []
trials = 1000
n = 6 #number of qubits in GHZ state
k = 2 #number of auxiliary qubits

for p in range(15):
    p = p/1000
    p_arr.append(p)
    all_sum_sc = 0
    for i in range(trials):
        all_sum_sc += sc_build_ghz_parallel(n-1,k,p)
    print("SC: " + str(all_sum_sc/trials))
    sc_arr.append(all_sum_sc/trials)
    
    all_sum_cn = 0
    for i in range(trials):
        all_sum_cn += cn_build_ghz_parallel(n,p)
    print("CN: " + str(all_sum_cn/trials))
    cn_arr.append(all_sum_cn/trials)

plt.show()
plt.scatter(p_arr,sc_arr,color = 'blue')
plt.scatter(p_arr,cn_arr, color = 'orange')
plt.xlabel("Error Probability")
plt.ylabel("Fidelity")
plt.legend(["Subgraph Complementation Algorithm","Factory node + Teleportation Algorithm"])
plt.axhline(y=0.75, color='red', linestyle='--')  # Draw a horizontal line at y=0.753

plt.show()