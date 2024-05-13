import networkx as nx
import numpy as np
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

#Subgraph Complementation Protocol
def sc(m,bit_str,err_str,p_e,last = False, reset = True): 
    n = len(bit_str)
    for i in range(1,n):
        cz(m,bit_str[0],bit_str[i],err_str,p_e)

    for i in range(1,n):
        lc(m,bit_str[i],err_str,p_e)

    for i in range(1,n):
        cz(m,bit_str[0],bit_str[i],err_str,p_e)

    lc(m,bit_str[0],err_str,p_e)

    if(last):
        for i in range(len(m)//2):
            z_meas(m,i,err_str,p_e)
        err_str = err_str[len(m)//2:]
        return int(err_str == [0]*(len(m)//2))

    if(reset):
        #Edge reset operation
        for i in range(1,n):
            cz(m,bit_str[0],bit_str[i],err_str,p_e)

        for i in range(1,n):
            lc(m,bit_str[i],err_str,p_e)

        for i in range(1,n):
            cz(m,bit_str[0],bit_str[i],err_str,p_e)

############################################################################################################

#Subgraph Complementation Algorithm (SC protocol) to create complete bipartite graph state
def sc_build_bip(n,p):
    m1 = [[0 for i in range(2*n)] for j in range(2*n)]
    err_str = [0 for i in m1]

    for i in range(n):
        cz(m1,i,n+i,err_str,p) #initialise graph state with bell pairs

    sc(m1,list(range(n)),err_str,p,last=False,reset=False) #list indices to do first SC on
    sc(m1,[0]+list(range(n))[n//2:],err_str,p,last=False,reset=False) #list indices to do second SC on
    return (sc(m1,list(range(n))[n//2:],err_str,p,last=True)) #list indices to do third SC on, thus creating a complete bipartite graph state

#Factory node + Teleportation Algorithm to create complete bipartite graph state
def cn_build_bip(n,p_e):
    m = [[0 for i in range(3*n)] for j in range(3*n)]
    err_str = [0 for i in m]
    
    for i in range(n):
        cz(m,i+n,i+2*n,err_str,p_e) #initialise graph state with bell pairs

    for i in range(n):
        for j in range(n//2,n):
            cz(m,i,j,err_str,p_e)
    
    for i in range(n):
        cz(m,i,i+n,err_str,p_e)

    for i in range(n):
        y_meas(m,i,err_str,p_e)
        y_meas(m,i+n,err_str,p_e)

    return 1 if err_str[2*n:] == [0]*n else 0


############################################################################################################

#Simulation of complete bipartite graph state generation using SC and CN protocols
# The following code generates graphs of error probability vs fidelity for SC and CN protocols, but this can be easily modified to generate graphs of number of qubits vs fidelity with small changes
sc_arr = []
cn_arr = []
p_arr = []
trials = 10000
n = 6
    
for p in range(0,30,2):
    p = p/10000
    p_arr.append(p)
    all_sum_sc = 0
    for i in range(trials):
        all_sum_sc += sc_build_bip(n,p)
    print("SC: " + str(all_sum_sc/trials))
    sc_arr.append(all_sum_sc/trials)
    
    all_sum_cn = 0
    for i in range(trials):
        all_sum_cn += cn_build_bip(n,p)
    print("CN: " + str(all_sum_cn/trials))
    cn_arr.append(all_sum_cn/trials)

plt.show()
plt.scatter(p_arr,sc_arr,color = 'blue')
plt.scatter(p_arr,cn_arr, color = 'orange')
plt.xlabel("Error probability")
plt.ylabel("Fidelity")
plt.legend(["Subgraph Complementation Algorithm","Factory node + Teleportation Algorithm"])
plt.axhline(y=0.75, color='red', linestyle='--')  # Draw a horizontal line at y=0.75

plt.show()
