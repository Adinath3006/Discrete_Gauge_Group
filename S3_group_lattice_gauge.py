import numpy as np
import random
from numba import jit,vectorize,njit
from numba.typed import listobject

# Matrix representations of S3 permutations
identity = np.matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

transposition_12 =np.matrix([[0, 1, 0],
                    [1, 0, 0],
                    [0, 0, 1]])

transposition_23 = np.matrix([[1, 0, 0],
                    [0, 0, 1],
                    [0, 1, 0]])

transposition_13 = np.matrix([[0, 0, 1],
                    [0, 1, 0],
                    [1, 0, 0]])

three_cycle_123 = np.matrix([[0, 1, 0],
                   [0, 0, 1],
                   [1, 0, 0]])

three_cycle_132 = np.matrix([[0, 0, 1],
                   [1, 0, 0],
                   [0, 1, 0]])

# Store the matrices in a list
s3_matrices = [identity, transposition_23, transposition_12, three_cycle_123, three_cycle_132, transposition_13]

basis = [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]

def basis_map(matrix):
    if (matrix == s3_matrices[0]).all():
        return np.transpose(np.matrix(basis[0]))
    elif (matrix == s3_matrices[1]).all():
        return np.transpose(np.matrix(basis[1]))
    elif (matrix == s3_matrices[2]).all():
        return np.transpose(np.matrix(basis[2]))
    elif (matrix == s3_matrices[3]).all():
        return np.transpose(np.matrix(basis[3]))
    elif (matrix == s3_matrices[4]).all():
        return np.transpose(np.matrix(basis[4]))
    else:
        return np.transpose(np.matrix(basis[5]))
    
def inverse_map(vector):
    if (vector == basis[0]).all():
        return 0
    elif (vector == basis[1]).all():
        return 1
    elif (vector == basis[2]).all():
        return 2
    elif (vector == basis[3]).all():
        return 3
    elif (vector == basis[4]).all():
        return 4
    else:
        return 5
  
def L_h_plus(h,basis):
    Lh = np.zeros([6,6])
    for i in range(6):
        temp = s3_matrices[h]*s3_matrices[i]
        bas = basis_map(temp)
        for j in range(6):
            bas_2 = np.matrix(basis[j])
            Lh[j][i] = bas_2*bas
    return Lh

def L_h_minus(h,basis):
    Lh = np.zeros([6,6])
    for i in range(6):
        temp = s3_matrices[i]*np.linalg.inv(s3_matrices[h])
        bas = basis_map(temp)
        for j in range(6):
            bas_2 = np.matrix(basis[j])
            Lh[j][i] = bas_2*bas
    return Lh

def display_matrix(h):
    #eigval, eigvec = np.linalg.eig(L_h_plus(h,basis))
    print(f'The L_h_plus matrix h = {h}')
    print(L_h_plus(h,basis))
    print('Corresponding eigenvectors')
    #print(eigvec)
         
def A0_g(g,state):
    for i in [3,0]:
        basis_i = np.transpose(np.matrix(basis[state[i]]))
        new_basis_i = L_h_plus(g,basis)*basis_i
        state[i] = inverse_map(np.transpose(new_basis_i))
    return state

def gauge_invariant_state(h):
    C = np.zeros([6,6,6,6])
    for i in range(6):
        for j in range(6):
            for k in range(6):
                for l in range(6):
                    state = [i,j,k,l]
                    new_state = A0_g(h,state)
                    if C[i][j][k][l] == 0:
                        C[i][j][k][l] = random.uniform(0,1)
                        C[new_state[0]][new_state[1]][new_state[2]][new_state[3]] = C[i][j][k][l]
    return C

mag_flux = []
for i in range(6):
    for j in range(6):
        for k in range(6):
            for l in range(6):
                matrix = s3_matrices[i]*s3_matrices[j]*s3_matrices[k]*s3_matrices[l]
                if (matrix == s3_matrices[0]).all():
                    mag_flux.append([i,j,k,l])
                    #print([i,j,k,l])

a = random.randint(0,5)
b = random.randint(0,5)
c = random.randint(0,5)
d = random.randint(0,5)
print(s3_matrices[a]*s3_matrices[b]*np.transpose(s3_matrices[c])*np.transpose(s3_matrices[d]))
print(np.trace(s3_matrices[a]*s3_matrices[b]*np.transpose(s3_matrices[c])*np.transpose(s3_matrices[d])))