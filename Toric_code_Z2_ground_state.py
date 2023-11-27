import numpy as np
import random

#------------------GAUGE GROUP---------------------

basis = [[1,0],[0,1]]

def basis_map(g):
    if g == 0:
        return np.transpose(np.matrix(basis[0]))
    else:
        return np.transpose(np.matrix(basis[1]))

def inverse_map(vector):
    if (vector == basis[0]).all():
        return 0
    else:
        return 1
    
def L_h_plus(h,basis):
    Lh = np.zeros([2,2])
    for i in range(2):
        temp = (h+i)%2
        bas = basis_map(temp)
        for j in range(2):
            bas_2 = np.matrix(basis[j])
            Lh[j][i] = bas_2*bas
    return Lh

#------------------LATTICE---------------------

# Number of links in the x direction with constant y is: n
n = 2
lattice = np.zeros([n,n],np.int32)

def v_links(v,lattice):
    m = lattice.shape[0]
    n = lattice.shape[1]
    i = v[0]
    j = v[1]
    if i == 0:
        if j == 0:
            links = [n*i+j,n*i+j+n-1,n**2+n*i+j,2*n**2-n]
        else:
            links = [n*i+j,n*i+j-1,n**2+n*i+j,2*n**2-n+j]
    else:
        if j == 0:
            links = [n*i+j,n*i+j+n-1,n**2+n*i+j,n**2+n*i-n+j]
        else:
            links = [n*i+j,n*i+j-1,n**2+n*i+j,n**2+n*i-n+j]
    return links
        
def A_v(v,state):
    sig_x = np.matrix([[0,1],[1,0]]) 
    links = v_links(v,lattice)
    for k in range(4):
        vec = sig_x*basis_map(state[links[k]])
        state[links[k]] = inverse_map(np.transpose(vec))
        
def B_p(p,state):
    l1 = p
    l2 = n**2+p
    l3 = (p+n)%(n**2)
    if (l2+1)%n == 0:
        l4 = l2-(p%n)
    else:
        l4 = l2-(p%n)+1
    sig_z = np.matrix([[1,0],[0,-1]]) 
    links = [l1,l2,l3,l4]
    for k in range(4):
        vec = sig_z*basis_map(state[links[k]])
        state[links[k]] = inverse_map(np.transpose(vec))
        
def dectobin(decimal,arr,dim,i):
    if(decimal > 0):
        arr[len(arr)-1-i] = decimal%2
        dectobin(int(decimal/2),arr,dim,i+1)
    return arr
        
def generate_basis_state(i,dim):
    arr = np.zeros(dim)
    return dectobin(i,arr,dim,0)
        
print(generate_basis_state(0,2*n**2))

def ground_state_coeff_A():
    c = np.zeros(2**n)
    

# Define random basis
state = []
for i in range(2*n**2):
    state.append(random.choice([0,1]))
    
state = [1, 1, 1, 1, 1, 1, 1, 1]
temp = state
print(state)

A_v([1,1],state)
print(v_links([1,1],lattice))
print(state)