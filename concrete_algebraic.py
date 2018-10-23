#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
 
from algebraic_structures import ComRing, PolyOverIntegralDomain as PolyID, \
    Field, PolyOverField, ComRingQuotientED, PolyOverGalois, EuclideanDomain
from math import gcd

from abc import abstractclassmethod

# =============================================================================
#https://stackoverflow.com/questions/4798654/modular-multiplicative-inverse-function-in-python
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
# =============================================================================

#TODO: Documentar

class Integer(EuclideanDomain):
    def __init__(self, int_ ):
        self.int_ = int_
    def div_mod(self, other):
        div_int_, mod_int_ = divmod(self.int_, other.int_)
        return Integer( div_int_ ), Integer( mod_int_ )
    @classmethod
    def euclidean_function(cls, obj):
        return abs( obj.int_ )
    
    def add(self, other):
        return Integer( self.int_ + other.int_ )
    def equals(self, other):
        return self.int_ == other.int_
    def mul(self, other):
        return Integer( self.int_ * other.int_ )
    
    def is_invertible(self):
        result = self.equals(self.__class__.one()) or self.equals(
                self.__class__.one().symmetric())
        return result
    
    @classmethod
    def one(cls):
        return Integer(1)
    def symmetric(self):
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


class Rational(Field):

    def __init__(self, p, q=1):
        if q == 0:
            raise ValueError('Denominator cannot be zero.')
        if p ==0:
            q = 1
        gcd_ = gcd(p, q)
        self.p = p // gcd_
        self.q = q // gcd_
    
    def add(self, other):
        #XXX: o denominador é self.q * other.q
        return Rational( self.p*other.q + self.q * other.p, self.q * other.q )
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
    def symmetric(self):
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
        if self.is_zero():
            raise ValueError('No inverse for zero.')
        return Rational(self.q, self.p)
    

class PolyOverQ( PolyOverField ):
    @classmethod
    def coefRing(cls):
        return Rational

class PolyOverZ(PolyID):
    
    @classmethod
    def coefRing(cls):
        return Integer
    
    @classmethod
    def from_int_coefs(cls, coefs):
        return PolyOverZ( [ Integer(x) for x in coefs ]  )
    
    
# =============================================================================

#TODO: docstring
class ZPrime(ComRingQuotientED, Field):    
    
    def copy(self):
        return self.__class__( self.rep )
    
    def inverse(self):
        if self.is_zero():
            raise ValueError('zero has no inverse')
        return self.__class__.from_int( 
                modinv( self.rep.int_, self.gen().int_ ) 
                )
    
    @classmethod
    def baseRing(cls):
        return Integer
    
    def represent(self):
        return str( self.rep )
    
    def __repr__(self):
        return self.represent()

    @classmethod
    def from_int(cls, i):
        return cls( Integer(i) )
        
def makeZp( p ):
    """
    Generate Zp from p prime.
    """
    #TODO: make sure p is prime
    class Zp( ZPrime ):
        @classmethod
        def gen(cls):
            return Integer( p )
    
    return Zp

def makeZpX( p ):
    Zp = makeZp( p )

    class ZpX(PolyOverGalois):
        @classmethod
        def coefRing(cls):
            return Zp
        @classmethod
        def from_int_coefs(cls, int_coefs):
            return cls( [ Zp.from_int(i) for i in int_coefs] )
        
    return ZpX
    
    
#TODO: PolyOverGalois quotient field which is a field

# =============================================================================

if __name__ == '__main__':
    Zp3 = makeZp( 3 )
    gal_z = Zp3.from_int( 5 )
    gal_w = Zp3.from_int( 2 )
    print( Zp3.can_reduce )
    print( gal_z, gal_w)
    print( gal_z == gal_w)

    Zp3X = makeZpX( 3 )
    f = Zp3X.from_int_coefs( [1, 50] )
    print(f)
    
    pol = Zp3X.from_dict( {0 : Zp3.from_int(5) , 5 : Zp3.from_int(1) } )
    print(pol)



    pol2 = Zp3X.from_dict( {1 : Zp3.from_int(5) , 2 : Zp3.from_int(1) } )
    other = pol + pol2

    print( gal_w / gal_w)





















