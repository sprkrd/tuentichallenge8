#!/usr/bin/env python3

from math import sqrt
from functools import reduce
from sympy.ntheory import factorint as factorize


class Lcg:
    def __init__(self, s, a, c, m):
        self.s, self.a, self.c, self.m = s, a, c, m

    def __call__(self, n):
        s, a, c, m = self.s, self.a, self.c, self.m
        ret = ( fastModPower(a, n, m)*s )%m + ( c*modSumOfPowers(a, n, m) )%m
        ret = ret%m
        return ret

    def next(self, x):
        a, c, m = self.a, self.c, self.m
        return ((a*x)%m + c)%m


def fastModPower(x, n, m):
    if n == 0:
        return 1
    mod_n_2 = fastModPower(x, n//2, m)
    ret = (mod_n_2*mod_n_2) % m
    if n%2:
        ret = (ret*x)%m
    return ret


def modSumOfPowers(a, n, m):
    """
    Computes 1 + a + a^2 + (...) + a^(n-1) (mod m)
    """
    return ( (fastModPower(a, n, m*(a-1))-1)//(a-1) ) % m


# base = [2,3,5,7]
# base_product = reduce(lambda acc,x: acc*x, base, 1)
# wheel = [i for i in range(base_product) if all(i%b != 0 for b in base)]
# def factorize(n):
    # factors = {}
    # for b in base:
        # while n%b == 0:
            # factors[b] = factors.get(b, 0) + 1
            # n = n // b
    # offset = 0
    # k = 1
    # d = wheel[k]
    # n_sqrt = int(sqrt(n))
    # while d <= n_sqrt:
        # while n%d == 0:
            # factors[d] = factors.get(d,0) + 1
            # n = n // d
            # n_sqrt = int(sqrt(n))
        # k += 1
        # if k == len(wheel):
            # k = 0
            # offset += base_product
        # d = offset + wheel[k]
    # if n > 1:
        # factors[n] = 1
    # return factors


def updateFactors(factors1, factors2):
    for p, k in factors2.items():
        factors1[p] = factors1.get(p, 0) + k


def multiplicativeOrder(a, factors_m):
    factors_totient = {}
    for p, k in factors_m.items():
        if k > 1:
            factors_totient[p] = factors_totient.get(p, 0) + k-1
        updateFactors(factors_totient, factorize(p-1))
    m = reduce(lambda acc,f: acc*f[0]**f[1], factors_m.items(), 1)
    order = reduce(lambda acc,f: acc*f[0]**f[1], factors_totient.items(), 1)
    # the order divides phi(m), so let's find smallest n | phi(m) s.t.
    # a^n = 1 (mod m)
    for p, k in factors_totient.items():
        for i in range(k):
            if fastModPower(a, order//p, m) == 1:
                order = order//p
    return order


def gcd(a,b):
    if b > a:
        a,b = b,a
    while b > 0:
        r = a%b
        a = b
        b = r
    return a


def firstRepetition(s,a,c,m):
    lcg = Lcg(s,a,c,m)
    x_0 = lcg(m) # make sure to skip aperiodic part
    m_ = m // gcd((a-1)*x_0+c, m)
    factors_m_a_1 = factorize(m_)
    updateFactors(factors_m_a_1, factorize(a-1))
    lam = multiplicativeOrder(a, factors_m_a_1)
    mu = 0
    while lcg(mu+lam) != lcg(mu):
        mu += 1
    return mu + lam + 1


def main():
    n = int(input())
    for case in range(1,n+1):
        s, a, c, m = map(int, input().split())
        print("Case #{}: {}".format(case, firstRepetition(s,a,c,m)))


if __name__ == "__main__":
    main()

