#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 22:33:52 2018

@author: henrique
"""

from algebraic_structures import ComRing, PolyComRing, Field, PolyOverField


class Integer(ComRing):
    def __init__(self, int_ ):
        self.int_ = int_
    
    def add(self, other):
        return Integer( self.int_ + other.int_ )
    def equals(self, other):
        return self.int_ == other.int_
    def mul(self, other):
        return Integer( self.int_ * other.int_ )
    @classmethod
    def one(cls):
        return Integer(1)
    def symetric(self):
        return Integer( -self.int_ )
    @classmethod
    def zero(cls):
        return Integer(0)
    def is_zero(self):
        return self.int_ == 0
    def is_one(self):
        return self.int_ == 1
    def copy(self):
        return Integer( self.int_.copy() )
    def represent(self):
        return str( self.int_ )
    def __repr__(self):
        return self.represent()


from math import gcd
class Rational(Field):

    def __init__(self, p, q):
        gcd_ = gcd(p, q)
        self.p = p // gcd_
        self.q = q // gcd_
    
    def add(self, other):
        return Rational( self.p*other.q + self.q * other.p, self.q * self.p )
    def equals(self, other):
        return self.p == other.p and self.q  == other.q
    def mul(self, other):
        return Rational( self.p * other.p, self.q * other.q )
    @classmethod
    def one(cls):
        return Rational(1, 1)
    
    @classmethod
    def zero(cls):
        return Rational(0, 1)    
    def symetric(self):
        return Rational( -self.p, self.q)
    
    def is_zero(self):
        return self.p == 0
    def is_one(self):
        return self.p == 1 and self.q == 1
    def copy(self):
        return Rational( self.p, self.q )
    def represent(self):
        return str( self.p ) + '/' + str( self.q )
    def __repr__(self):
        return self.represent()
    def inverse(self):
        return Rational(self.q, self.p)
    

class PolyOverQ( PolyOverField ):
    @classmethod
    def coefRing(cls):
        return Rational

class PolyOverZ(PolyComRing):
    
    @classmethod
    def coefRing(cls):
        return Integer
    
    @classmethod
    def from_int_coefs(cls, coefs):
        return PolyOverZ( [ Integer(x) for x in coefs ]  )


f = PolyOverZ.from_int_coefs( [1, 2, 1] )
g = PolyOverZ.from_int_coefs( [1, -2, 3] )

print('f = ', f)
print('g = ', g)


h = f + g
p = f*g
print('f + g = h =', h)
print('f * g = p = ', p)

# =============================================================================

print('\nOver Q\n')
phi = PolyOverQ( [Rational(1,2), Rational(2,3)] )
print( phi )
phi2 = phi**2

