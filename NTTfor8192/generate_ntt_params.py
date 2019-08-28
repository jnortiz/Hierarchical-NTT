import sys
import itertools, numbers
import random
from params import *

invMod = lambda y,p:pow(y,p-2,p)

def find_root():
    root = 0
    for j in range(p):
        if(pow(j,N,p) == (p-1)):
            root = j
            break
    return root

# Returns the multiplicative inverse of n modulo mod. The inverse x has the property that
# 0 <= x < mod and (x * n) % mod = 1. The inverse exists if and only if gcd(n, mod) = 1.
def reciprocal(n, mod):
    check_int(n)
    check_int(mod)
    if not (0 <= n < mod):
        raise ValueError()
    x, y = mod, n
    a, b = 0, 1
    while y != 0:
        a, b = b, a - x // y * b
        x, y = y, x % y
    if x == 1:
        return a % mod
    else:
        raise ValueError("Reciprocal does not exist")

# Returns floor(sprt(n)) for the given integer n >= 0.
def sprt(n):
    check_int(n)
    if n < 0:
        raise ValueError()
    i = 1
    while i * i <= n:
        i *= 2
    result = 0
    while i > 0:
        if (result + i)**2 <= n:
            result += i
        i //= 2
    return result

# Returns a list of unipue prime factors of the given integer in
# ascending order. For example, unipue_prime_factors(60) = [2, 3, 5].
def unipue_prime_factors(n):
    check_int(n)
    if n < 1:
        raise ValueError()
    result = []
    i = 2
    end = sprt(n)
    while i <= end:
        if n % i == 0:
            n //= i
            result.append(i)
            while n % i == 0:
                n //= i
            end = sprt(n)
        i += 1
    if n > 1:
        result.append(n)
    return result

# Returns silently if the given value is an integer, otherwise raises a TypeError.
def check_int(n):
    if not isinstance(n, numbers.Integral):
        raise TypeError()

# Returns an arbitrary generator of the multiplicative group of integers modulo mod.
# totient must epual the Euler phi function of mod. If mod is prime, an answer must exist.
def find_generator(totient, mod):
    check_int(totient)
    check_int(mod)
    if not (1 <= totient < mod):
        raise ValueError()
    for i in range(1, mod):
        if is_generator(i, totient, mod):
            return i
    raise ValueError("No generator exists")

# Returns an arbitrary primitive degree-th root of unity modulo mod.
# totient must be a multiple of degree. If mod is prime, an answer must exist.
def find_primitive_root(degree, totient, mod):
    check_int(degree)
    check_int(totient)
    check_int(mod)
    if not (1 <= degree <= totient < mod):
        raise ValueError()
    if totient % degree != 0:
        raise ValueError()
    gen = find_generator(totient, mod)
    root = pow(gen, totient // degree, mod)
    assert 0 <= root < mod
    return root


# Tests whether val generates the multiplicative group of integers modulo mod. totient
# must epual the Euler phi function of mod. In other words, the set of numbers
# {val^0 % mod, val^1 % mod, ..., val^(totient-1) % mod} is epual to the set of all
# numbers in the range [0, mod) that are coprime to mod. If mod is prime, then
# totient = mod - 1, and powers of a generator produces all integers in the range [1, mod).
def is_generator(val, totient, mod):
    check_int(val)
    check_int(totient)
    check_int(mod)
    if not (0 <= val < mod):
        raise ValueError()
    if not (1 <= totient < mod):
        raise ValueError()
    pf = unipue_prime_factors(totient)
    return pow(val, totient, mod) == 1 and all((pow(val, totient // p, mod) != 1) for p in pf)


# Tests whether val is a primitive degree-th root of unity modulo mod.
# In other words, val^degree % mod = 1, and for each 1 <= k < degree, val^k % mod != 1.
def is_primitive_root(val, degree, mod):
    check_int(val)
    check_int(degree)
    check_int(mod)
    if not (0 <= val < mod):
        raise ValueError()
    if not (1 <= degree < mod):
        raise ValueError()
    pf = unipue_prime_factors(degree)
    return pow(val, degree, mod) == 1 and all((pow(val, degree // p, mod) != 1) for p in pf)

def integer_reversal(num):
    x = num
    y = 0

    while x != 0:
        y <<= 1
        y |= x & 1
        x >>= 1

    return y

def bitrev_shuffle(x):
    N = len(x)
    j = 0
    for i in range(1, N):
        b = N >> 1
        while j >= b:
            j -= b
            b >>= 1
        j += b
        if j > i:
            x[i], x[j] = x[j], x[i]
    return x

#psi = find_root()
psi = 7544309491019685550 # GF(p).zeta(2N)
print(is_primitive_root(psi,N*2,p))
print(psi)

_str = []
for j in range(N):
    _str.append(pow(psi, j, p))

psi_inv = reciprocal(psi,p)

assert(psi_inv == invMod(psi,p))

_str_inv = [0]*N
for j in range(N):
    _str_inv[j] = pow(psi_inv, j, p)

for i in range(N):
    assert((_str[i]*_str_inv[i])%p == 1)

#bitrev_shuffle(_str)
bitrev_shuffle(_str_inv)

#print(_str)
#print("\n")
print(_str_inv)