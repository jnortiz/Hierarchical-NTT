invMod = lambda y,p:pow(y,p-2,p)

# p = 9223372036801560577
# Nr = 64
# Nc = 128
# N = Nr*Nc

# g = 8492750271230630006
# g_inv = invMod(g,p)

# psi = 1679062545781875027
# psi_inv = invMod(psi, p)

# root_Nr = 617500568322588441
# root_Nr_inv = invMod(root_Nr, p)

# root_Nc = 8769891463669861743
# root_Nc_inv = invMod(root_Nc, p)

p = 12289
primitive_root_of_p = 11
Nr = 2
Nc = 4
N = Nr*Nc # N = 8

g = pow(primitive_root_of_p, (p-1)//N, p) #primitive 8-th root of unity in Z_p
g_inv = invMod(g,p)

psi = pow(primitive_root_of_p, (p-1)//(2*N), p) #primitive 16-th root of unity in Z_p
psi_inv = invMod(psi, p)

root_Nr = pow(primitive_root_of_p, (p-1)//Nr, p) #primitive 2-th root of unity
root_Nr_inv = invMod(root_Nr, p)

root_Nc = pow(primitive_root_of_p, (p-1)//Nc, p) #primitive 4-th root of unity
root_Nc_inv = invMod(root_Nc, p)

assert(pow(g, N, p) == 1)
assert(pow(g_inv, N, p) == 1)

assert(pow(psi, 2*N, p) == 1)
assert(pow(psi_inv, 2*N, p) == 1)