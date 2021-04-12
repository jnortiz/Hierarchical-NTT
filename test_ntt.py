#!/usr/env/python3
#coding: utf-8

import random
import unittest
from params import *
from ntt import ntt_gs, intt_ct, poly_mul, rewritten_ntt_gs, rewritten_intt_ct, poly_mul_rewritten


N_RUNS = 1000


def schoolbook_mul(a, b):
    assert len(a) == len(b)

    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            v = a[i] * b[j] * (-1) ** (int((i + j) // float(N)))
            c[(i + j) % N] = (c[(i + j) % N] + v) % p
    
    return c


def gen_polynomial_modp(length):
    x = []
    for i in range(length):
        x.append(random.randrange(0, p))
    return x


class TestNTT(unittest.TestCase):

    def test_transformation(self):

        print("\n NTT transformation")

        for i in range(N_RUNS):
            a = gen_polynomial_modp(N)
            b = intt_ct(ntt_gs(a))
            self.assertEqual(a, b)

    def test_poly_mul(self):

        print("\n Polynomial multiplication using NTT-GS and INTT-CT")

        for i in range(N_RUNS):
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)
            c_prime = schoolbook_mul(a, b)
            c = poly_mul(a, b)
            self.assertEqual(c_prime, c)

    def test_transformation_rewritten(self):

        print("\n Rewritten NTT transformation")

        for i in range(N_RUNS):
            a = gen_polynomial_modp(N)
            b = rewritten_intt_ct(rewritten_ntt_gs(a))
            self.assertEqual(a, b)

    def test_poly_mul_rewritten(self):

        print("\n Polynomial multiplication using rewritten NTT-GS and INTT-CT")

        for i in range(N_RUNS):
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)
            c_prime = schoolbook_mul(a, b)
            c = poly_mul_rewritten(a, b)
            self.assertEqual(c_prime, c)


if __name__ == '__main__':
    unittest.main()
