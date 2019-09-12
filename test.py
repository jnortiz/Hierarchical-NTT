import unittest
from params import *
from ntt import *

N_RUNS = 10**1

def remove_zeros(y):
    
    x = list(y)

    while(x[0] == 0):
        x.remove(0)

    while(x[len(x)-1] == 0):
        x.pop()

    return x

class TestNTT(unittest.TestCase):

    def test_transformation(self):
        print("\nHierarchical NTT transformation")
        
        for i in range(N_RUNS):
            
            a = gen_polynomial_modp(N)

            b = hierarchical_intt(hierarchical_ntt(a))

            self.assertEqual(a, b)

    def test_mul(self):        
        print("\nHierarchical polynomial multiplication")

        for _ in range(1):

            a = [1]*N
            b = [1]*N
            
            c = hierarchical_ntt_mul(a, b)

            print(c)
            #self.assertEqual(ab,c)

if __name__ == '__main__':
    unittest.main()