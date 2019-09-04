#!/usr/env/python3
#coding: utf-8
import sys
import random
from numpy import reshape, array, polymul
from math import log,floor
from params import *

invMod = lambda y,p:pow(y,p-2,p)
bitreversal = lambda n,width:int('{:0{width}b}'.format(n, width=width)[::-1], 2)

def schoolbook_mul(a, b):
    assert len(a) == len(b)
    
    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            v = a[i]*b[j]*(-1)**(int((i+j)//float(N)))
            c[(i+j) % N] = (c[(i+j) % N] + v) % p
    return c   

def transpose(x):
    return array(x).reshape(Na, Nb).transpose().reshape(Na * Nb).tolist()

def gen_polynomial_modp(length):
    x = []
    for i in range(length):
        x.append(random.randrange(0,p))
    return x

def hierarchical_ntt(x):
    assert(len(x) == N)

    s = list(x)

    s = [y*pow(psi,i,p)%p for i,y in enumerate(s)]  # Negative wrapped convolution      
    s = transpose(s)
    s = sum([list(ntt_simple(s[i*Na:i*Na+Na],Na,root_Na)) for i in range(Nb)], [])
    s = [y*pow(g,(i // Na) * (i % Na),p)%p for i,y in enumerate(s)]    
    s = sum([list(ntt_simple(s[i::Na],Nb,root_Na)) for i in range(Nb)],[])

    return s

def hierarchical_intt(x):
    assert(len(x) == N)

    s = list(x)

    s = sum([list(intt_simple(s[i*Nb:i*Nb+Nb],Nb,root_Na_inv)) for i in range(Na)],[])
    s = [y*pow(g_inv,(i // Na) * (i % Na),p)%p for i,y in enumerate(s)]    
    s = sum([list(intt_simple(s[i::Na],Nb,root_Na_inv)) for i in range(Nb)], [])    
    s = transpose(s)    
    s = [y*pow(psi_inv,i,p)%p for i,y in enumerate(s)] # Negative wrapped convolution   

    return s  

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
    b_ntt = hierarchical_ntt(b)

    c_ntt = poly_mul_pointwise(a_ntt, b_ntt)

    c = hierarchical_intt(c_ntt)

    return c