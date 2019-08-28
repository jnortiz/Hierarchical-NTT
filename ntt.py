#!/usr/env/python3
#coding: utf-8
import sys
import unittest
import time
import random
from numpy import transpose
from math import log,floor
from params_ntt import *

p = 9223372036801560577
N = 8192

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
    Na = 64
    Nb = 128
    assert(Na*Nb == N)
    assert(len(poly_x) == N)

    x = poly_x
    x_ntt = sum([ntt128(x[i*Nb:i*Nb+Nb]) for i in range(Na)], [])

    print(len(x_ntt))

    for i in range(Na):
        for j in range(Nb):
            x_ntt[i*Nb+j] = x_ntt[i*Nb+j]*pow(factor,i*Nb+j,p)

    x_ntt2 = transpose([ntt64(x_ntt[i::Nb]) for i in range(Na)]).reshape(Na*Nb).tolist()[0]

    return x_ntt2

def hierarchical_intt(poly_x_ntt): # The matrix is transposed!
    Na = 64
    Nb = 128
    assert(Na*Nb == N)
    assert(len(poly_x_ntt) == N)

    x_ntt = poly_x_ntt
    x_ntt2 = sum([intt64(x_ntt[i*Na:i*Na+Na]) for i in range(Nb)], [])

    for i in range(Nb):
        for j in range(Na):
            x_ntt2[i*Na+j] = x_ntt2[i*Na+j]*pow(factor_inv,i*Na+j,p)

    x = transpose([intt128(x_ntt2[i::Na]) for i in range(Nb)]).reshape(Nb*Na).tolist()[0]

    return x

def poly_mul_pointwise(x_ntt, y_ntt):
    xy = [x*y for x,y in zip(x_ntt, y_ntt)]
    return xy

def ntt64(x):
    N = 64    
    a = x
    t = N
    m = 1
    while m < N:
        t >>= 1
        for i in range(m):
            j1 = 2*i*t
            j2 = j1+t-1
            S = pow(root_64,bitreversal(m+i,6),p)
            for j in range(j1,j2+1):
                U = a[j]
                V = a[j+t]*S
                a[j] = (U+V)%p
                a[j+t] = (U-V)%p
        m<<=1
    return a

def intt64(x):
    N = 64
    a = list(x)
    t = 1
    m = N
    while m > 1:
        j1 = 0
        h = m >> 1
        for i in range(h):
            j2 = j1 + t - 1
            S = pow(rootinv_64,bitreversal(h+i,6),p)
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

def ntt128(x):
    N = 128
    a = list(x)
    t = N
    m = 1
    while m < N:
        t >>= 1
        for i in range(m):
            j1 = 2*i*t
            j2 = j1+t-1
            S = pow(root_128,bitreversal(m+i,7),p)
            for j in range(j1,j2+1):
                U = a[j]
                V = a[j+t]*S
                a[j] = (U+V)%p
                a[j+t] = (U-V)%p
        m<<=1
    return a

def intt128(x):
    N = 128
    a = list(x)
    t = 1
    m = N
    while m > 1:
        j1 = 0
        h = m >> 1
        for i in range(h):
            j2 = j1 + t - 1
            S = pow(rootinv_128,bitreversal(h+i,7),p)
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