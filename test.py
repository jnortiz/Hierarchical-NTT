import unittest
from params import *
from ntt import *

N_RUNS = 1

def schoolbook_mul(a, b):
    assert len(a) == len(b)
    
    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            v = a[i]*b[j]*(-1)**(int((i+j)//float(N)))
            c[(i+j) % N] = (c[(i+j) % N] + v) % p
    return c    

class TestNTT(unittest.TestCase):

    def test_transformation(self):
        print("\nHierarchical NTT transformation")
        
        for i in range(N_RUNS):
            a = gen_polynomial_modp(N)

            b = hierarchical_intt(hierarchical_ntt(a))

            self.assertEqual(a, b)

    def test_mul(self):        
        print("\nHierarchical polynomial multiplication")

        for _ in range(N_RUNS):
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)

            ab = schoolbook_mul(a,b)
            c = hierarchical_ntt_mul(a, b)

            self.assertEqual(ab,c)

if __name__ == '__main__':
    unittest.main()