import unittest
from params import *
from ntt import *

N_RUNS = 10**2

def remove_zeros(y):
    
    x = list(y)

    while(x[0] == 0):
        x.remove(0)

    while(x[len(x)-1] == 0):
        x.pop()

    return x

class TestNTT(unittest.TestCase):

    def test_ntt_simple(self):
        print("\nStraightforward NTT")

        for _ in range(N_RUNS):

            a = gen_polynomial_modp(N)        
            b = intt_simple(ntt_simple(a,N,g),N,g_inv)
            
            self.assertEqual(a,b)

    def test_negative_wrapped_convolution(self):
        print("\nPolynomial multiplication using straightforward NTT with negative wrapped convolution")

        for _ in range(N_RUNS):

            a = gen_polynomial_modp(N)        
            b = gen_polynomial_modp(N)

            a_hat = [x*pow(psi,i,p)%p for i,x in enumerate(a)]
            b_hat = [x*pow(psi,i,p)%p for i,x in enumerate(b)]

            a_ntt = ntt_simple(a_hat,N,g)
            b_ntt = ntt_simple(b_hat,N,g)

            c_ntt = [(x*y) for x,y in zip(a_ntt, b_ntt)]

            c_hat = intt_simple(c_ntt,N,g_inv)

            c = [x*pow(psi_inv,i,p)%p for i,x in enumerate(c_hat)]

            ab = schoolbook_mul(a,b)

            self.assertEqual(ab,c)

    def test_ntt_simple_mul(self):
        print("\nPolynomial multiplication using straightforward NTT without polynomial reduction")

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

    def test_transformation(self):
        print("\nHierarchical NTT transformation")
        
        for _ in range(N_RUNS):        
            a = gen_polynomial_modp(N)
            b = hierarchical_intt(hierarchical_ntt(a))
            self.assertEqual(a, b)

    def test_mul(self):        
        print("\nPolynomial multiplication using <<hierarchical>> straightforward NTT")

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

    def test_mul_negacyclic(self):        
        print("\nPolynomial multiplication using <<hierarchical>> straightforward NTT with negacyclic")

        for _ in range(N_RUNS):
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)

            ab = schoolbook_mul(a,b)            
            
            c = hierarchical_ntt_mul(a, b)

            self.assertEqual(ab,c)            

if __name__ == '__main__':
    unittest.main()