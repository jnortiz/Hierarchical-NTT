from numpy import reshape, array
from math import log

from ntt import rewritten_ntt_gs as ntt, rewritten_intt_ct as intt
from params import *


def hierarchical_ntt(x):

    assert(len(x) == N)

    s = list(x)

    # Forward negative wrapped convolution
    # psi is a power of the primitive root
    s = [s[i] * pow(psi, i, p) % p for i in range(N)]

    # Forward hierarchical NTT
    s = sum([list(ntt(s[i:: Nc])) for i in range(Nc)], [])
    s = [s[i] * pow(g, (i // Nr) * bit_reverse(i %
                    Nr, int(log(Nr, 2))), p) % p for i in range(N)]
    s = sum([list(ntt(s[i:: Nr])) for i in range(Nr)], [])

    return s


def hierarchical_intt(x):

    assert(len(x) == N)

    s = list(x)

    # Inverse hierarchical NTT
    s = sum([list(intt(s[i * Nc: i * Nc + Nc])) for i in range(Nr)], [])
    s = [s[i] * pow(g_inv, bit_reverse(i // Nc, int(log(Nr, 2)))
                    * (i % Nc), p) % p for i in range(N)]
    s = sum([list(intt(s[i:: Nc])) for i in range(Nc)], [])

    # This transposition format the output as row-major
    s = array(s).reshape(Nc, Nr).transpose().reshape(Nr * Nc).tolist()

    # Inverse negative wrapped convolution
    s = [(s[i] * pow(psi_inv, i, p)) % p for i in range(N)]

    return s


def poly_mul_pointwise(x_ntt, y_ntt):
    return [(x * y) % p for x, y in zip(x_ntt, y_ntt)]


def hierarchical_ntt_mul(a, b):

    a_ntt = hierarchical_ntt(a)
    b_ntt = hierarchical_ntt(b)

    c_ntt = poly_mul_pointwise(a_ntt, b_ntt)

    c = hierarchical_intt(c_ntt)

    return c
