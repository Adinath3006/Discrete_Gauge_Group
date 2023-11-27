import numpy as np
import math
import itertools

#------------------GAUGE GROUP---------------------

def generate_permutation_matrix(perm, n):
    matrix = np.zeros((n, n))
    for i in range(n):
        matrix[i, perm[i]-1] = 1
    return matrix

def generate_symmetric_group_representations(n):
    representations = []
    elements = list(range(1, n+1))
    permutations = list(itertools.permutations(elements))
    
    for perm in permutations:
        matrix = generate_permutation_matrix(perm, n)
        representations.append(np.matrix(matrix))
    
    return representations

N = int(input("Enter the value of N: "))
representations = generate_symmetric_group_representations(N)

# Print and store the representations
def display_representation():
    for i, representation in enumerate(representations):
        print(f"Representation {i+1}:")
        print(representation,'\n')

# Create a basis
def generate_basis(n):
    basis = []
    for i in range(math.factorial(n)):
        basis_element = np.zeros(math.factorial(n))
        basis_element[i] = 1
        basis.append(basis_element)
    return basis

# Create a mapping between the elements of the S_n group and the basis
def generate_basis_map(matrix,n):
    for i in range(math.factorial(n)):
        if (matrix == generate_symmetric_group_representations(n)[i]).all():
            return np.transpose(np.matrix(generate_basis(n)[i]))

# Create the left action operator   
def L_h_plus(n,h):
    basis = generate_basis(n)
    m = math.factorial(n)
    sn_matrices = generate_symmetric_group_representations(n)
    Lh = np.zeros([m,m])
    for i in range(m):
        temp = np.matrix(sn_matrices[h])*np.matrix(sn_matrices[i])
        bas = generate_basis_map(temp,n)
        for j in range(m):
            bas_2 = np.matrix(basis[j])
            Lh[j][i] = bas_2*bas
    return Lh

# Create the right action operator
def L_h_minus(n,h):
    basis = generate_basis(n)
    m = math.factorial(n)
    sn_matrices = generate_symmetric_group_representations(n)
    Lh = np.zeros([m,m])
    for i in range(m):
        temp = np.matrix(sn_matrices[i])*np.linalg.inv(np.matrix(sn_matrices[h]))
        bas = generate_basis_map(temp,n)
        for j in range(m):
            bas_2 = np.matrix(basis[j])
            Lh[j][i] = bas_2*bas
    return Lh

#------------------LATTICE---------------------

# Number of links in the x direction with constant y is: n
n = 2
lattice = np.zeros([n,n],np.int32)

# Gives the set of links emanating from the vertex i
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

# Gives the set of links on the plaquette p
def p_plaquette(p,n):
    l1 = p
    l2 = n**2+p
    l3 = (p+n)%(n**2)
    if (l2+1)%n == 0:
        l4 = l2-(p%n)
    else:
        l4 = l2-(p%n)+1
    links = [l1,l4,l3,l2]
    return links

#------------------COMPUTE GROUND STATE---------------------

def ground_state(n,N):
    gnd_states = []
    dim_state = 2*n**2
    states = list(itertools.product(range(math.factorial(N)), repeat=dim_state))
    m = len(states)
    for i in range(m):
        check = 1
        for p in range(n**2):
            link = p_plaquette(p,n)
            prod = representations[states[i][link[0]]]*representations[states[i][link[1]]]*np.linalg.inv(representations[states[i][link[2]]])*np.linalg.inv(representations[states[i][link[3]]])
            if (prod == np.matrix(representations[0])).all(): 
                check *= 1
            else:
                check *= 0
        if check == 1:
            print(states[i])
            gnd_states.append(states[i])
    return gnd_states

ground_state(n,N)