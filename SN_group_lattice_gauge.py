import numpy as np
import math
import itertools

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

n = int(input("Enter the value of n: "))
representations = generate_symmetric_group_representations(n)

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

# Display the operators
def display_matrix(n,h):
    #eigval, eigvec = np.linalg.eig(L_h_plus(n,h))
    print(f'The L_h_plus matrix h = {h}')
    print(L_h_plus(n,h))
    print('Corresponding eigenvectors')
    #print(eigvec)
    

mag_flux = []
for i in range(math.factorial(n)):
    for j in range(math.factorial(n)):
        for k in range(math.factorial(n)):
            for l in range(math.factorial(n)):
                matrix = representations[i]*representations[j]*representations[k]*representations[l]
                if (matrix == np.matrix(representations[0])).all():
                    mag_flux.append([i,j,k,l])
                    print([i,j,k,l])
                    
print(len(mag_flux))