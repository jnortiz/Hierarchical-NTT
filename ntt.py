#!/usr/env/python3
#coding: utf-8
import sys
import unittest
import time
import random
from numpy import transpose
from math import log,floor
from params_ntt import *

#python -m unittest ntt.TestNTT.test_transformation

p = 9223372036801560577
N = 8192
Na = 64
Nb = 128

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

def gen_polynomial_modp(length):
    x = []
    for i in range(length):
        x.append(random.randrange(0,p))
    return x

def hierarchical_ntt(poly_x):
    assert(Na*Nb == N)
    assert(len(poly_x) == N)

    x = poly_x
    x_ntt = sum([ntt(x[i*Nb:i*Nb+Nb],128,root_128) for i in range(Na)], [])

    for i in range(Na):
        for j in range(Nb):
            x_ntt[i*Nb+j] = x_ntt[i*Nb+j]*pow(factor,i*Nb+j,p)

    x_ntt2 = transpose([ntt(x_ntt[i::Nb],64,root_64) for i in range(Nb)]).reshape(Na*Nb).tolist()

    return x_ntt2

def hierarchical_intt(poly_x_ntt): # The matrix is transposed!
    assert(Na*Nb == N)
    assert(len(poly_x_ntt) == N)

    x_ntt = poly_x_ntt
    x_ntt2 = sum([intt(x_ntt[i*Na:i*Na+Na],64,rootinv_64) for i in range(Nb)], [])

    for i in range(Nb):
        for j in range(Na):
            x_ntt2[i*Na+j] = x_ntt2[i*Na+j]*pow(factor_inv,i*Na+j,p)

    x = transpose([intt(x_ntt2[i::Na],128,rootinv_128) for i in range(Na)]).reshape(Nb*Na).tolist()

    return x

def poly_mul_pointwise(x_ntt, y_ntt):
    xy = [x*y for x,y in zip(x_ntt, y_ntt)]
    return xy

def ntt(x,N,root):
    a = x
    log_N = int(log(N,2))
    t = N
    m = 1
    while m < N:
        t >>= 1
        for i in range(m):
            j1 = 2*i*t
            j2 = j1+t-1
            S = pow(root,bitreversal(m+i,log_N),p)
            for j in range(j1,j2+1):
                U = a[j]
                V = a[j+t]*S
                a[j] = (U+V)%p
                a[j+t] = (U-V)%p
        m<<=1
    return a

def intt(x,N,rootinv):    
    a = list(x)
    log_N = int(log(N,2))
    t = 1
    m = N
    while m > 1:
        j1 = 0
        h = m >> 1
        for i in range(h):
            j2 = j1 + t - 1
            S = pow(rootinv,bitreversal(h+i,log_N),p)
            for j in range(j1, j2+1):
                U = a[j]
                V = a[j+t]
                a[j] = (U+V)%p
                a[j+t] = ((U-V)*S)%p
            j1 = j1 + 2*t
        t <<= 1
        m >>= 1

    for j in range(N):
        a[j] = (a[j]*invMod(N,p))%p
    return a

def hierarchical_ntt_mul(a, b):
    
    a_ntt = hierarchical_ntt(a)
    b_ntt = hierarchical_ntt(b)

    c_ntt = poly_mul_pointwise(a_ntt, b_ntt)

    c = hierarchical_intt(c_ntt)

    return c 

class TestNTT(unittest.TestCase):

    def test_transformation(self):
        print("\nTesting NTT Gentleman-Sande")
        a = gen_polynomial_modp(N)

        start_time = time.time()        
        b = hierarchical_intt(hierarchical_ntt(a))
        end_time = time.time()

        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(a, b)

    def test_mul(self):        
        print("\nPolynomial multiplication using NTT Gentleman-Sande")

        a = gen_polynomial_modp(N)
        b = gen_polynomial_modp(N)

        start_time = time.time()        
        c = hierarchical_ntt_mul(a, b)
        end_time = time.time()
      
        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(schoolbook_mul(a,b), c)

if __name__ == '__main__':
    unittest.main()