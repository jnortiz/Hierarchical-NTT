invMod = lambda y,p:pow(y,p-2,p)

p = 9223372036801560577
Nr = 64
Nc = 128
N = Nr*Nc

g = 8492750271230630006
g_inv = invMod(g,p)

psi = 1679062545781875027
psi_inv = invMod(psi, p)

root_Nr = 617500568322588441
root_Nr_inv = invMod(root_Nr, p)

root_Nc = 8769891463669861743
root_Nc_inv = invMod(root_Nc, p)