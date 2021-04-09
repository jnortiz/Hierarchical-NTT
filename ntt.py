from math import log

from params import *


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


def poly_mul_pointwise(x_ntt, y_ntt):
    return [(x * y) % p for x, y in zip(x_ntt, y_ntt)]


# gentleman-sande ntt
def ntt_gs(x):

    r = list(x)
    N_local = len(r)
    g_local = pow(primitive_root_of_p, (p - 1) // N_local, p)

    m = N_local // 2
    while m >= 1:
        for j in range(m):
            a = pow(g_local, (j * N_local) // (2 * m), p)
            i = j
            while i < N_local:
                u = r[i]
                v = r[i + m]
                r[i] = (u + v) % p
                r[i + m] = a * (u - v) % p
                i = i + 2 * m
        m >>= 1

    return bitrev_shuffle(r)


# cooley-tukey intt
def intt_ct(x):

    r = list(x)
    r = bitrev_shuffle(r)

    N_local = len(r)
    g_inv_local = invMod(pow(primitive_root_of_p, (p - 1) // N_local, p), p)

    m = 1
    while m < N_local:
        for j in range(m):
            a = pow(g_inv_local, (j * N_local) // (2 * m), p)
            i = j
            while i < N_local:
                u = r[i]
                v = a * r[i + m] % p
                r[i] = (u + v) % p
                r[i + m] = (u - v) % p
                i = i + 2 * m
        m <<= 1

    return [(x * invMod(N_local, p)) % p for x in r]


# gentleman-sande ntt
def rewritten_ntt_gs(x):

    r = list(x)
    N_local = len(r)
    g_local = pow(primitive_root_of_p, (p - 1) // N_local, p)

    for s in range(int(log(N_local, 2))):
        m = N_local // (2 << s)
        for l in range(N_local//2):
            j = (2 * m * l) // N_local
            i = j + (l % (N_local // (2 * m))) * (2 * m)
            a = pow(g_local, j * (N_local >> (int(log(N_local, 2)) - s)), p)
            u = r[i]
            v = r[i + m]
            r[i] = (u + v) % p
            r[i + m] = a * (u - v) % p

    return r


# cooley-tukey intt
def rewritten_intt_ct(x):

    r = list(x)

    N_local = len(r)
    g_local = pow(primitive_root_of_p, (p - 1) // N_local, p)

    m = 1
    for s in range(int(log(N_local, 2))):
        for l in range(N_local // 2):
            j = (2 * m * l) // N_local
            i = j + (l % (N_local // (2 * m))) * (2 * m)
            a = pow(g_local, (N_local - j) * (N_local >> (s + 1)), p)
            u = r[i]
            v = r[i + m]
            r[i] = (u + a * v) % p
            r[i + m] = (u - a * v) % p
        m <<= 1

    return [(x * invMod(N_local, p)) % p for x in r]


def poly_mul(a, b):

    a = [ai * pow(psi, i, p) % p for i, ai in enumerate(a)]
    b = [bi * pow(psi, i, p) % p for i, bi in enumerate(b)]

    a_ntt = ntt_gs(a)
    b_ntt = ntt_gs(b)

    c_ntt = poly_mul_pointwise(a_ntt, b_ntt)

    c = intt_ct(c_ntt)

    c = [ci * pow(psi_inv, i, p) % p for i, ci in enumerate(c)]

    return c


def poly_mul_rewritten(a, b):

    a = [ai * pow(psi, i, p) % p for i, ai in enumerate(a)]
    b = [bi * pow(psi, i, p) % p for i, bi in enumerate(b)]

    a_ntt = rewritten_ntt_gs(a)
    b_ntt = rewritten_ntt_gs(b)

    c_ntt = poly_mul_pointwise(a_ntt, b_ntt)

    c = rewritten_intt_ct(c_ntt)

    c = [ci * pow(psi_inv, i, p) % p for i, ci in enumerate(c)]

    return c
