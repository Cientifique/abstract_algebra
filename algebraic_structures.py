#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

# =============================================================================

from abc import ABC, abstractmethod, abstractclassmethod
 
# =============================================================================

def rep_default(x):
    type_ = type(x)
    module = type_.__module__
    qualname = type_.__qualname__
    return f"<{module}.{qualname} object at {hex(id(x))}>"

# =============================================================================

class ComRing(ABC):
    """
    Commutative Ring R(+,*)
    """
    
    @abstractmethod
    def __init__(self, *args, **kwargs):
        #Constructor
        pass
    
    @abstractclassmethod
    def zero(cls):
        """
        Creates the zero element.
        """
        pass
    
    @abstractclassmethod
    def one(cls):
        """
        Creates the one element.
        """
        pass
    
    @abstractmethod
    def add(self, other):
        """
        Returns self + other.
        """
        pass
    
    @abstractmethod
    def mul(self, other):
        """
        Returns self * other.
        """
        pass

    @abstractmethod
    def copy(self):
        """
        Returns a copy.
        """
        pass
    
    def natural_power(self, k):
        """
        Returns self ^k for a non-negative integer k.
        """
        
        #Makes sure it is an integer
        k = int(k)
        if k < 0:
            raise ValueError('k Must non-negative integer')
        if k == 0:
            return ComRing.one()
        elif k == 1:
            return self.copy()
        else:
            result = self
            for k_ in range(k):
                result = result.mul(result)
            return result

    def power(self, k):
        """
        This methood can be made more general.
        """
        return self.natural_power(k)
    
    @abstractmethod
    def symetric(self):
        """
        Returns symetric, so that self + symetric == 0
        """
        pass
    
    @abstractmethod
    def is_zero(self):
        """
        Returns boolean to indicate if self is the zero of the ring.
        """
        pass
    
    @abstractmethod
    def is_one(self):
        """
        Returns boolean to indicate if self is the zero of the ring.
        """
        pass
    
    @abstractmethod
    def equals(self, other):
        """
        Returns boolean indicating if self == other.
        """
        pass
    
    #OPERATOR OVERRIDE
    
    def __eq__(self, other):
        return self.equals(other)
    
    def __add__(self, other):
        return self.add(other)
    
    def __neg__(self):
        return self.symetric()
    
    def __sub__(self, other):
        return self.add( other.symetric() )
    
    def __mul__(self, other):
        return self.mul( other )
    
    def __pow__(self, k):
        return self.power(k)
            
    def __radd__(self, other):
        return other.add(self)
    
# =============================================================================

class EuclideanDomain(ComRing):
    """
    Euclidean Domain.
    """
    @abstractclassmethod
    def euclidean_function(cls, x):
        """
        Euclidean function f so that in the Euclidean Domain R:
            
            If a and b are in R and b is nonzero, then there are q and r in R 
            such that a = bq + r and either r = 0 or f(r) < f(b).
        
        """
        pass
    
    @abstractmethod
    def div_mod(self, other):
        """
        Division with remainder. Returns quotient, remainder.
        """
        pass
    
    def quotient(self, other):
        """
        Returns the quotient of the euclidean division in Zp[X].
        """
        result, _ = self.div_mod(other)
        return result

    def mod(self, other):
        """
        Returns the remainder of the euclidean division.
        """
        _, result = self.div_mod(other)
        return result
    
    #OPERATOR OVERRIDE
    
    def __floordiv__(self, other):
        """
        Overriding. Allows to write f // g.
        """
        result = self.quotient(other)
        return result        

    def __mod__(self, other):
        """
        Overriding. Allows to write f % g.
        """
        return self.mod(other)

# =============================================================================

class Field(EuclideanDomain):
    
    @abstractmethod
    def inverse(self):
        """
        Returns multiplicative inverse.
        """
        pass
    
    def div(self, other):
        """
        Divides self by other: self * other.inverse()
        """
        return self.mul( other.inverse() )
    
    #Now we can use negative powers
    def power(self, k):
        k = int(k)
        if k >= 0:
            return self.natural_power(k)
        else:
            return self.inverse().natural_power(k)
    
    @classmethod
    def euclidean_function(cls, x):
        return cls.one()
    
    def div_mod(self, other):
        return self.div(other), self.zero()
    
    def __truediv__(self, other):
        return self.div(other)

# =============================================================================

class PolyComRing(ComRing):
    
    #Assuming ring of coefs is an integral domain, maybe with class attribute
    
    @abstractclassmethod
    def coefRing(cls):
        pass
    
    @staticmethod
    def _remove_trailing_zeros(coefs):
        k = 0
        for k, value in enumerate( coefs[::-1] ):
            if not value.is_zero():
                break
        lst_no_trailing_zeroes = coefs if k == 0 else coefs[:-k]
        return lst_no_trailing_zeroes
    
    @classmethod
    def _validate_coefs(cls, coefs):
        """
        Makes sure coefs are valid.
        """
        if any( not isinstance(x , cls.coefRing() ) for x in coefs):
            raise ValueError('All coefficients must belong to coefRing')
    
    def __init__(self, coefs, validate_coefs = False, remove_trailing_zeroes = False):
        if validate_coefs:
            PolyComRing.check_coefs( coefs )
        if remove_trailing_zeroes:
            coefs_ = PolyComRing._remove_trailing_zeros( coefs )
        else:
            coefs_ = coefs
        self.coefs = coefs_
    
    def is_zero(self):
        return self.coefs == []
    
    def is_one(self):
        if self.degree() != 1:
            return False
        else:
            return self.constant_coef().is_one()
    
    def degree(self):
        if self.is_zero():
            return -1
        else:
            return len(self.coefs) - 1
    
    def leading_coef(self):
        return self.coefs[ self.degree() - 1]
    
    def constant_coef(self):
        if self.degree() == -1:
            return 0
        else:
            return self.coefs[0]
    
    def evaluate(self, x0):
        """
        Evaluates the polynomial at x0.
        """
        return sum(  
                (
                        coef * x0**k 
                        for coef, k in enumerate( self.coefs )
                        if not x0.is_zero()
                        ),
                self.coefRing().zero()
                )
        
    def remove(self, idx):
        """
        Removes coef from poly.
        """
        if idx == self.degree():
            self.coefs = self.coefs[:-1]
        else:
            self.coefs[idx] = self.coefRing().zero()
    
    @classmethod
    def zero(cls):
        return cls([])
    
    @classmethod
    def one(cls):
        return cls( [ cls.coef_one() ] )
    
    def symetric(self):
        new = self.__class__( [ -coef for coef in self.coefs] )
        return new
    
    def equals(self, other):
        return self.coefs == other.coefs
    
    def copy(self):
        copy_coefs = self.coefs.copy()
        copy = self.__class__( copy_coefs )    
        return copy
    
    def add(self, other):
        common_coefs = [ coef1 + coef2 for coef1, coef2 in zip(self.coefs, other.coefs) ]
        if self.degree() == other.degree():
            res_coefs = self._remove_trailing_zeros( common_coefs )
        else:
            if self.degree() >= other.degree():
                big_ = self
            else:
                big_ = other
                
            res_coefs = common_coefs + big_.coefs[ len(common_coefs) : ]
        result = self.__class__( res_coefs )
        return result
        
    def _mul_monomial(self, mono_coef, mono_power):
        res_coefs = [ self.coefRing().zero() for x in range( mono_power)]
        for coef in self.coefs:
            res_coefs.append( coef * mono_coef )
        result_ = self.__class__( res_coefs )
        return result_
    
    def mul(self, other):
       
        result = sum(
                ( 
                        self._mul_monomial(mono_coef, mono_power)
                        for mono_power, mono_coef in enumerate( other.coefs )
                ),
                self.zero()
                )
        return result
    
    
    def __repr__(self):
        
        """
        Returns a string with the usual representation of the polynomial.
        It is in LaTeX compatible form.
        """
        if hasattr( self.coefRing(), 'represent'):
            #Trivial case for the zero polynomial
            if self.is_zero():
                return '0'
            else:
                deg = self.degree()
                #Degree zero just returns the constant coeficient
                if deg == 0:
                    return str( self.coefs[0] )
                
                else:
                    poly_str_lst = []
                    #Special formatting for constant coef and x coef
                    if not self.coefs[0].is_zero():
                        poly_str_lst.append( str(self.coefs[0]) )
                    first_coef = self.coefs[1]
                    if not first_coef.is_zero():
                        first_coef_ = '' if first_coef.is_one() else first_coef
                        poly_str_lst.append( '%sx' % first_coef_ )
                    if deg == 1:
                        poly_str = ' + '.join( poly_str_lst )
                    else:
                        #adds powers to the rest, and hides coef if it is 1
                        for power, coef in zip( range(2, deg+1), self.coefs[2:] ):
                            if not coef.is_zero():
                                coef_ = '' if coef.is_one() else coef
                                poly_str_lst.append( '%sx^{%s}' % ( coef_, power) )
                        poly_str = ' + '.join( poly_str_lst )
                    return poly_str
        else:
            return rep_default(self)

# =============================================================================
      
class PolyOverField(PolyComRing, EuclideanDomain):
    """
    Ring of Polynomials over Field.
    """
    
    def euclidean_function(cls, g ):
        return g.degree()
        
    #TODO
    def div_mod(self, other):
        pass
    


# =============================================================================


#TODO: Function that generates the class Zp(Field) for given p prime
#TODO: Function that generates the class ZpX(Field) for given p prime (use above)
#TODO: Function that takes Ring R and a method is_in_I and creates quotient ring Q = R/I
#TODO: Class PolyOverGalois(PolyOverField) with is_irreducible implemented



        
        
        
        
        


