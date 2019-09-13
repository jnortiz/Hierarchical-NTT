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
Nr = 2
Nc = 4
N = Nr*Nc # N = 8

g = 8246 # 8-th root of unity in Z_p
g_inv = invMod(g,p)

psi = 4134 #16-th root of unity in Z_p
psi_inv = invMod(psi, p)

root_Nr = 12288 #2-th root of unity
root_Nr_inv = invMod(root_Nr, p)

root_Nc = 1479 #4-th root of unity
root_Nc_inv = invMod(root_Nc, p)