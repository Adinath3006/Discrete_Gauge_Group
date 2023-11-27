import numpy as np
import random

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

def display_matrix(h):
    eigval, eigvec = np.linalg.eig(L_h_plus(h,basis))
    print(f'The L_h_plus matrix h = {h}')
    print(L_h_plus(h,basis))
    print('Corresponding eigenvectors')
    print(eigvec)

def A0_g(g,state):
    for i in range(4):
        basis_i = np.transpose(np.matrix(basis[state[i]]))
        new_basis_i = L_h_plus(g,basis)*basis_i
        state[i] = inverse_map(np.transpose(new_basis_i))
    return state

def prntC(C):
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    print(C[i][j][k][l])

def gauge_invariant_state(h):
    C = np.zeros([2,2,2,2])
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    state = [i,j,k,l]
                    new_state = A0_g(h,state)
                    if C[i][j][k][l] == 0:
                        C[i][j][k][l] = random.uniform(0,1)
                        C[new_state[0]][new_state[1]][new_state[2]][new_state[3]] = C[i][j][k][l]
    return C
                        
def print_state(h):
    coeff = gauge_invariant_state(h)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    if i == 1 and j == 1 and k == 1 and l == 1:
                        print(f'{coeff[i][j][k][l]}|{i}{j}{k}{l}>')
                    else:
                        print(f'{coeff[i][j][k][l]}|{i}{j}{k}{l}> + ',end = '')

#print_state(0)
display_matrix(1)