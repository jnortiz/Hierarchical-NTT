from params import *


def bit_reverse(n, width): return int(
    '{:0{width}b}'.format(n, width=width)[::-1], 2)


def log_2(n): return n.bit_length()


def transpose(x, Nr, Nc):
    return array(x).reshape(Nr, Nc).transpose().reshape(Nr*Nc).tolist()


def bitrev_shuffle(x):
    N = len(x)
    j = 0
    for i in range(1, N):
        b = N >> 1
        while j >= b:
            j -= b
            b >>= 1
        j += b
        if j > i:
            x[i], x[j] = x[j], x[i]
    return x


def ntt_gs(x):

    r = list(x)
    N_local = len(r)
    g_local = pow(primitive_root_of_p, (p-1)//N_local, p)

    m = N_local//2
    while m >= 1:
        for j in range(m):
            a = pow(g_local, (j*N_local)//(2*m), p)
            i = j
            while i < N_local:
                u = r[i]
                v = r[i+m]
                r[i] = (u + v) % p
                r[i+m] = a*(u - v) % p
                i = i + 2 * m
        m >>= 1

    return bitrev_shuffle(r)


def intt_ct(x):

    r = list(x)
    r = bitrev_shuffle(r)

    N_local = len(r)
    g_inv_local = pow(pow(primitive_root_of_p, (p-1)//N_local, p), p-2, p)

    m = 1
    while m < N_local:
        for j in range(m):
            a = pow(g_inv_local, (j*N_local)//(2*m), p)
            i = j
            while i < N_local:
                u = r[i]
                v = a * r[i+m] % p
                r[i] = (u + v) % p
                r[i + m] = (u - v) % p
                i = i + 2*m
        m <<= 1

    return [(x * pow(N_local, p-2, p)) % p for x in r]


def poly_mul(a, b):

    a = [ai * pow(psi, i, p) % p for i, ai in enumerate(a)]
    b = [bi * pow(psi, i, p) % p for i, bi in enumerate(b)]

    a_ntt = ntt_gs(a)
    b_ntt = ntt_gs(b)

    for i in range(N):
        c_ntt = [(x*y) % p for x, y in zip(a_ntt, b_ntt)]

    c = intt_ct(c_ntt)

    c = [ci * pow(psi_inv, i, p) % p for i, ci in enumerate(c)]

    return c
