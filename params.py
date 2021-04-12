def bit_reverse(n, width): return int(
    '{:0{width}b}'.format(n, width=width)[::-1], 2)


def invMod(y, p): return pow(y, p - 2, p)


# Real parameters
p = 0xFFFFFFFF00000001
primitive_root_of_p = 7
Nr = 4
Nc = 8
N = Nr * Nc

# Toy example using NewHope parameters:

# p = 12289
# primitive_root_of_p = 11
# Nr = 32
# Nc = 32
# N = Nr*Nc

g = pow(primitive_root_of_p, (p - 1) // N, p)
g_inv = invMod(g, p)

psi = pow(primitive_root_of_p, (p - 1) // (2 * N), p)
psi_inv = invMod(psi, p)

assert(pow(g, N, p) == 1)
assert(pow(g_inv, N, p) == 1)

assert(pow(psi, 2 * N, p) == 1)
assert(pow(psi_inv, 2 * N, p) == 1)
