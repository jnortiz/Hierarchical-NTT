from ntt import ntt_gs, intt_ct
from params import *
from numpy import reshape, array

def hierarchical_ntt(x):

    assert(len(x) == N)

    s = list(x)

    # Forward negative wrapped convolution
    s = [s[i]*pow(psi, i, p) % p for i in range(N)]

    # This transposition format the input as column-major (as in Fortran language, which was used by David H. Bailey)
    s = array(s).reshape(Nr, Nc).transpose().reshape(Nr*Nc).tolist()

    # Forward hierarchical NTT
    s = sum([list(ntt_gs(s[i*Nr:i*Nr+Nr]))
            for i in range(Nc)], [])
    s = [s[i]*pow(g, (i//Nr)*(i % Nr), p) % p for i in range(N)]
    s = array(s).reshape(Nc, Nr).transpose().reshape(Nr*Nc).tolist()
    s = sum([list(ntt_gs(s[i*Nc:i*Nc+Nc]))
            for i in range(Nr)], [])

    return s


def hierarchical_intt(x):

    assert(len(x) == N)

    s = list(x)

    # Inverse hierarchical NTT
    s = sum([list(intt_ct(s[i*Nc:i*Nc+Nc]))
            for i in range(Nr)], [])
    s = [s[i]*pow(g_inv, (i//Nc)*(i % Nc), p) % p for i in range(N)]
    s = array(s).reshape(Nr, Nc).transpose().reshape(Nr*Nc).tolist()
    s = sum([list(intt_ct(s[i*Nr:i*Nr+Nr]))
            for i in range(Nc)], [])

    # This transposition format the output as row-major
    s = array(s).reshape(Nc, Nr).transpose().reshape(Nr*Nc).tolist()

    # Inverse negative wrapped convolution
    s = [s[i]*pow(psi_inv, i, p) % p for i in range(N)]

    return s


def poly_mul_pointwise(x_ntt, y_ntt):
    return [(x*y) % p for x, y in zip(x_ntt, y_ntt)]


def hierarchical_ntt_mul(a, b):
    
    a_ntt = hierarchical_ntt(a)
    b_ntt = hierarchical_ntt(b)
    
    c_ntt = poly_mul_pointwise(a_ntt, b_ntt)
    
    c = hierarchical_intt(c_ntt)

    return c
