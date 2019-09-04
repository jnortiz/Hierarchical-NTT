#!/usr/env/python3
#coding: utf-8
import sys
import unittest
import time
import random
from numpy import reshape, array, polymul
from math import log,floor

N_RUNS = 10**3
p = 12289
g = 7311 # 64-th root of unity in Z_p
Na = 8
Nb = 8
N = Na*Nb

invMod = lambda y,p:pow(y,p-2,p)

rootNa = 8246 #8-th root of unity
rootNainv = invMod(rootNa,p)

bitreversal = lambda n,width:int('{:0{width}b}'.format(n, width=width)[::-1], 2)

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

    s = transpose(s)
    s = sum([ntt_simple(s[i*Na:i*Na+Na],Na,rootNa) for i in range(Nb)], [])
    s = [y*pow(g,(i // Na) * (i % Na),p)%p for i,y in enumerate(s)]
    s = sum([ntt_simple(s[i::Na],Nb,rootNa) for i in range(Nb)],[])

    return s

def hierarchical_intt(x):
    assert(len(x) == N)

    s = list(x)

    s = sum([intt_simple(s[i*Nb:i*Nb+Nb],Nb,rootNainv) for i in range(Na)],[])
    s = [y*pow(invMod(g,p),(i // Na) * (i % Na),p)%p for i,y in enumerate(s)]
    s = (sum([intt_simple(s[i::Na],Nb,rootNainv) for i in range(Nb)], []))
    s = transpose(s)

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

    s = []
    for k in range(N):
        s.append(
            invMod(N,p)*sum(
                [x[j]*pow(rootinv, j*k, p) % p for j in range(N)]
                ) % p
            )    

    return s

def ntt(x,N,root):
    a = list(x)
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

def remove_zeros(y):
    
    x = list(y)

    while(x[0] == 0):
        x.remove(0)

    while(x[len(x)-1] == 0):
        x.pop()

    return x

class TestNTT(unittest.TestCase):

    def test_ntt_simple(self):
        print("\nStraighforward NTT")

        g = 7311 # N-th root of unity in Z_p
        g_inv = invMod(g,p)

        for _ in range(N_RUNS):        
            a = gen_polynomial_modp(N)        
            b = intt_simple(ntt_simple(a,N,g),N,g_inv)
            self.assertEqual(a,b)

    def test_ntt_simple_mul(self):
        print("\nPolynomial multiplication using straighforward NTT without polynomial reduction")

        g = 7311 # N-th root of unity in Z_p
        g_inv = invMod(g,p)

        for _ in range(N_RUNS):        
            a = gen_polynomial_modp(N//2)        
            b = gen_polynomial_modp(N//2)

            ab = list(polymul(a,b)%p)

            a = a + [0]*(N//2)
            b = b + [0]*(N//2)

            a_ntt = ntt_simple(a,N,g)
            b_ntt = ntt_simple(b,N,g)

            c_ntt = [(x*y) for x,y in zip(a_ntt, b_ntt)]

            c = intt_simple(c_ntt,N,g_inv)

            c = remove_zeros(c)
            ab = remove_zeros(ab)

            self.assertEqual(ab,c)

    def test_ntt(self):
        print("\nPolynomial multiplication using NTT Gentleman-Sande with negacyclic convolution")

        g = 12149 # 2N-th root of unity in Z_p
        g_inv = invMod(g,p)

        for _ in range(N_RUNS):        
            a = gen_polynomial_modp(N//2)        
            b = gen_polynomial_modp(N//2)

            ab = list(polymul(a,b)%p)

            a = a + [0]*(N//2)
            b = b + [0]*(N//2)

            a_ntt = ntt(a,N,g)
            b_ntt = ntt(b,N,g)

            c_ntt = [(x*y) for x,y in zip(a_ntt, b_ntt)]

            c = intt(c_ntt,N,g_inv)

            c = remove_zeros(c)
            ab = remove_zeros(ab)

            self.assertEqual(ab,c)

    def test_transformation(self):
        print("\nHierarchical NTT transformation")
        
        for _ in range(N_RUNS):        
            a = gen_polynomial_modp(N)
            b = hierarchical_intt(hierarchical_ntt(a))
            self.assertEqual(a, b)

    def test_mul(self):        
        print("\nPolynomial multiplication using <<hierarchical>> straighforward NTT")

        for _ in range(N_RUNS):        
            a = gen_polynomial_modp(N//2)
            b = gen_polynomial_modp(N//2)
            
            ab = list(polymul(a,b)%p)

            a = a + [0]*(N//2)
            b = b + [0]*(N//2)

            c = hierarchical_ntt_mul(a, b)

            c = remove_zeros(c)
            ab = remove_zeros(ab)

            self.assertEqual(ab,c)

if __name__ == '__main__':
    unittest.main()