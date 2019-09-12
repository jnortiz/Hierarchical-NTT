#!/usr/env/python3
#coding: utf-8
import sys
import random
from numpy import reshape, array, polymul
from math import log,floor
from params import *

def transpose(x, Nr, Nc):
    return array(x).reshape(Nr, Nc).transpose().reshape(Nr*Nc).tolist()

def gen_polynomial_modp(length):
    x = []
    for i in range(length):
        x.append(random.randrange(0,p))
    return x

def hierarchical_ntt(x):
    
    assert(len(x) == N)

    rows = Nr
    cols = Nc

    s = list(x)
    s = [s[i]*pow(psi,i,p)%p for i in range(N)]  # Negative wrapped convolution
    s = transpose(s, rows, cols)

    rows = Nc
    cols = Nr
  
    s = sum([list(ntt_simple(s[i*cols:i*cols+cols],cols,root_Nr)) for i in range(rows)], [])
    s = [s[i]*pow(g,(i//cols)*(i%cols),p)%p for i in range(N)]
    s = transpose(s, rows, cols)

    rows = Nr
    cols = Nc

    s = sum([list(ntt_simple(s[i*cols:i*cols+cols],cols,root_Nc)) for i in range(rows)], [])
     
    return s

def hierarchical_intt(x):

    assert(len(x) == N)

    rows = Nr
    cols = Nc

    s_ntt = list(x)
    s_ntt = sum([list(intt_simple(s_ntt[i*cols:i*cols+cols],cols,root_Nc_inv)) for i in range(rows)], [])
    s_ntt = [s_ntt[i]*pow(g_inv,(i//cols)*(i%cols),p)%p for i in range(N)]
    s_ntt = transpose(s_ntt, rows, cols)

    rows = Nc
    cols = Nr

    s_ntt = sum([list(intt_simple(s_ntt[i*cols:i*cols+cols],cols,root_Nr_inv)) for i in range(rows)], [])
    s_ntt = transpose(s_ntt, rows, cols)
    s_ntt = [s_ntt[i]*pow(psi_inv,i,p)%p for i in range(N)]  # Negative wrapped convolution  

    return s_ntt  

def poly_mul_pointwise(x_ntt, y_ntt):
    xy = [(x*y)%p for x,y in zip(x_ntt, y_ntt)]
    return xy

def ntt_simple(x,N,root):
    assert(pow(root,N,p) == 1) #root should be the n-th root in unity in Z_q
    assert(len(x) == N)

    s = []
    for k in range(N):
        s.append(sum([x[j]*pow(root, j*k, p) for j in range(N)]) % p)    

    return s

def intt_simple(x,N,rootinv):
    assert(pow(rootinv,N,p) == 1) #root should be the n-th root in unity in Z_q

    s = []
    for k in range(N):
        s.append(
            invMod(N,p)*sum(
                [x[j]*pow(rootinv, j*k, p) % p for j in range(N)]
                ) % p
            )    

    return s

def hierarchical_ntt_mul(a, b):
    
    a_ntt = hierarchical_ntt(a)
    assert(hierarchical_intt(a_ntt) == a)
    
    b_ntt = hierarchical_ntt(b)
    assert(hierarchical_intt(b_ntt) == b)

    c_ntt = poly_mul_pointwise(a_ntt, b_ntt)

    c = hierarchical_intt(c_ntt)

    return c