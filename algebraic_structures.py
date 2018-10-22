#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

# =============================================================================

from abc import ABC, abstractmethod, abstractclassmethod
 
# =============================================================================

def rep_default(x):
    """
    Returns the default representation of a python object
    "<{module}.{qualname} object at {hex(id(x))}>"
    
    example <function f at 0x7faeac0211e0> 
    """
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
            raise ValueError('k must be a non-negative integer')
        if k == 0:
            return self.one()
        elif k == 1:
            return self.copy()
        else:
            result = self
            for k_ in range(k):
                result = result.mul( result )
            return result

    def power(self, k):
        """
        
        Arguments:
            k - non-negative integer
        
        Returns:
            self raised to the power k.
        """
        return self.natural_power(k)
    
    @abstractmethod
    def symmetric(self):
        """
        Returns symmetric, so that self + symmetric == 0
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
        """ Allows to write self == other """
        return self.equals(other)
    
    def __add__(self, other):
        """ Allows to write self + other """
        return self.add(other)
    
    def __neg__(self):
        """ Allows to write -self """
        return self.symmetric()
    
    def __sub__(self, other):
        """ Allows to write self - other """
        return self.add( other.symmetric() )
    
    def __mul__(self, other):
        """ Allows to write self * other """
        return self.mul( other )
    
    def __pow__(self, k):
        """ Allows to write self ** k """
        return self.power(k)
            
    def __radd__(self, other):
        """ Allows to use the sum builtin """
        return other.add( self )
    
# =============================================================================

class EuclideanDomain(ComRing):
    """
    Euclidean Domain.
    """
    
    @abstractclassmethod
    def euclidean_function(cls, x):
        """
        Euclidean function f : R\{0} ---> N in the Euclidean Domain R that
        satisfies:
            
            - if a and b are in R and b is nonzero then there are unique q and r in R 
                such that a = bq + r and either r = 0 or f(r) < f(b);
            - if a and b are nonzero then f(a) <= f(ab)
        
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
        Returns the quotient of the euclidean division.
        """
        result, _ = self.div_mod(other)
        return result

    def mod(self, other):
        """
        Returns the remainder of the euclidean division.
        """
        _, result = self.div_mod(other)
        return result
    
    def isdivisible(self, other):
        """
        Returns bool - other divides self
        """
        return ( self % other ).is_zero()
    
    #OPERATOR OVERRIDE
    
    def __floordiv__(self, other):
        """
        Allows to write self // other.
        """
        result = self.quotient(other)
        return result        

    def __mod__(self, other):
        """
        Allows to write self % other.
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
        """
        Raise self to integer power. For negative powers it multiplies
        by the inverse of self raised to -k ( > 0).
        
        Arguments:
           k - int 
         
        Returns:
            self raised to the power k.
        """
        k = int(k)
        if k >= 0:
            return self.natural_power(k)
        else:
            return self.inverse().natural_power( -k )
    
    @classmethod
    def euclidean_function(cls, x):
        return 1
    
    #Division with remainder is just normal division in a Field
    def div_mod(self, other):
        return self.div(other), self.zero()
    
    def __truediv__(self, other):
        """ Allows to write self / other """
        return self.div(other)

# =============================================================================

class PolyOverIntegralDomain(ComRing):
    
    """
    Polynomial Ring over Commutative Ring which is an Integral Domain.
    
    Arguments:
        coefs - List of elements of the coefRing starting with constant
            term up to highest degree. Zero polynomial has coefs = []
        validate_coefs - boolean default False. If True makes sure all coefs
            belong to coefRing.
        remove_trailing_zeroes - boolean default False. If True removes all
            trailing zeroes from the coefs.
    
    """
    
    @abstractclassmethod
    def coefRing(cls):
        """
        Returns the class that corresponds to the coeficients of the Polynomial.
        """
        pass
    
    @staticmethod
    def _remove_trailing_zeros(coefs):
        """
        Removes trailling zeroes from the list of coefficients.
        
        Arguments:
            coefs - List of elements of the coefRing.
        """
        
        k = 0
        for k, value in enumerate( coefs[::-1] ):
            if not value.is_zero():
                break
        
        res_ = coefs if k == 0 else coefs[:-k]
        
        if len(res_) == 1:
            result = [] if res_[0].is_zero() else res_
        else:
            result = res_
        
        return result
    
    @classmethod
    def _validate_coefs(cls, coefs):
        """
        Validates if all coefficients belong to coefRing.
        """
        if any( not isinstance(x , cls.coefRing() ) for x in coefs):
            raise ValueError('All coefficients must belong to coefRing')
    
    @classmethod
    def from_dict(cls, dct, validate_coefs = False):
        """
        In case we want to instantiate a polynomial using a dict.
        Makes it easier to create polynomials such as 1 + x^100.
        It is assumed the dict has the form {..., power : coefficient, ...}, 
        where power is int and coefficient is from coefRing.
        Coefficients are assumed to be zero if they have no corresponding key.
        
        Arguments:
            dct - Dict of form { power : coef }
            validate_coefs - boolean like constructor.
        """
        if len(dct) == 0:
            result = cls([])
        else:
            if validate_coefs:
                cls._validate_coefs( dict.values() )
            
            powers = [ int(x) for x in dct.keys() ]
            degree = max(powers)
            
            # An m-degree polynomial has m+1 coefficients
            coefs = [ cls.coefRing().zero() for x in range(degree+1) ]
            
            # fills in the coefficients
            for i in powers:
                coefs[i] = dct[ i ]
            result = cls(coefs)
        
        return result
    
    @classmethod
    def zero(cls):
        return cls( [] )
    
    @classmethod
    def one(cls):
        return cls( [ cls.coefRing().one() ] )
    
    #Constructor is concrete
    def __init__(self, coefs, 
                 validate_coefs = False, remove_trailing_zeroes = False):
        
        if validate_coefs:
            self._validate_coefs( coefs )
        if remove_trailing_zeroes:
            coefs_ = self._remove_trailing_zeros( coefs )
        else:
            coefs_ = coefs
            
        #Keeps the coefs of the polynomial
        self.coefs = coefs_
        
    def typesetter(func):
        """
        Decorator.
        Allows binary operations to accept several input types: class instances 
        of polynomials (the default assumption), elements of the ring of 
        coefficients, lists and dicts.
        """
        def new_func(self, other):
            
            # in this case there is nothing to do, the default methods already 
            #assumes that the input is a class instance
            if isinstance(other, self.__class__):
                result = func(self, other)
            
            # instantiates a constant polynomial if other is an element of coefRing
            elif isinstance(other, self.__class__.coefRing()):
                coefs = [other]
                other = self.__class__(coefs, validate_coefs=False, 
                                       remove_trailing_zeroes=True)
                result = func(self, other)
            
            # instantiates a polynomial from the list
            elif type(other) == list:
                other = self.__class__(other, validate_coefs=True, 
                                       remove_trailing_zeroes=True)
                result = func(self, other)
            
            # instantiates a polynomial from the dict
            elif type(other) == dict:
                other = self.__class__.from_dict(other, validate_coefs = False)
                result = func(self, other)
            
            else:
                error_msg = 'Not a valid polynomial for binary operation.'
                raise TypeError(error_msg)
            return result
        
        return new_func
            
    
    def is_zero(self):
        return self.coefs == []
    
    #TODO: testar
    def is_one(self):
        if self.degree() != 0:
            return False
        else:
            return self.equals(self.__class__.one())
    
    def degree(self):
        """
        Returns the degree of the polynomial. If it is the zero polynomial
            then it returns -1.
        """
        if self.is_zero():
            return -1
        else:
            return len(self.coefs) - 1
    
    def leading_coef(self):
        """ Returns the leading coefficient of the polynomial"""
        if self.is_zero():
            return self.coefRing().zero()
        else: 
            return self.coefs[ self.degree() ]
    
    def constant_coef(self):
        if self.is_zero():
            return self.coefRing().zero()
        else:
            return self.coefs[0]
    
    def evaluate(self, x0):
        """
        Evaluates the polynomial at x0.
        
        Arguments:
            x0 - coefRing object
        Returns
            coefRing object
        """
        return sum( (
                coef * x0**k
                for coef, k in enumerate( self.coefs ) if not x0.is_zero()
                    ),
                self.coefRing().zero()
                )
        
    def remove(self, idx):
        """
        Removes coef from polynomial making it zero.
        
        Arguments:
            idx - int order of the coef. If < 0 makes no changes.
        """
        
        if idx >= 0 and not self.is_zero():
            if idx == self.degree():
                self.coefs = self.coefs[:-1]
            else:
                self.coefs[idx] = self.coefRing().zero()
    
    def symmetric(self):
        new = self.__class__( [ -coef for coef in self.coefs] )
        return new
    
    @typesetter
    def equals(self, other):
        return self.coefs == other.coefs
    
    def copy(self):
        copy_coefs = self.coefs.copy()
        copy = self.__class__( copy_coefs )    
        return copy
    
    @typesetter
    def add(self, other):
        common_coefs = [ 
                coef1 + coef2 for coef1, coef2 in zip(self.coefs, other.coefs)
                ]
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
        """
        Multiplies polynomial by monomial.
        
        Arguments:
            mono_coef - coefficient of monomial.
            mono_power - power of x for coefficient
        """
        
        #If any of them is zero return zero
        if mono_coef.is_zero() or self.is_zero():
            return self.zero()
        
        #The first mono_power coefficients are zero
        res_coefs = [ self.coefRing().zero() for x in range( mono_power)]
        #Then makes the product coef by coef
        for coef in self.coefs:
            res_coefs.append( coef * mono_coef )
        #Creates new object
        result_ = self.__class__( res_coefs )
        return result_
    
    @typesetter
    def mul(self, other):
       
        result = sum( ( 
                        self._mul_monomial(mono_coef, mono_power)
                        for mono_power, mono_coef in enumerate( other.coefs )
                ),
                self.zero()
                )
        return result
    
    
    def __repr__(self):
        
        """
        Returns a string with the usual representation of the polynomial 
            if the coefRing has a represent methood. Otherwise returns
            the default python representation of an object.
        It is in LaTeX compatible form for integers, reals, rationals...
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
      
class PolyOverField(PolyOverIntegralDomain, EuclideanDomain):
    """
    Ring of Polynomials over Field. It is an EuclideanDomain.
    """
    
    def euclidean_function(cls, g ):
        """
        Returns the degree of the polynomial.
        """
        return g.degree()
    
    #TODO: typesetter, talvez definir abstract method na classe PolyOverIntegralDomain
    # e pô-lo lá
    def div_mod(self, other):
        """
        Polynomial euclidean division.
        
        Pseudocode:
        
        Begin
            q := 0
            r := a
            d := deg(b)
            c := lc(b)
            while deg(r) ≥ d do
                s := lc(r)/c xdeg(r)−d
                q := q + s
                r := r − sb
            end do
            return (q, r)
        end.

        """
        #Cannot divide by zero
        if other.is_zero():
            raise ValueError('other must be non-zero')
        
        #Initialize quotient and remainder
        q = self.zero()
        # At each step self = d × q + r
        r = self
        
        #degree of divisor
        d = other.degree()
        
        #leading coef of divisor
        c = other.leading_coef()
        
        while r.degree() >= d:
            
            # Divide the leading terms
            s_coef = r.leading_coef() / c # times x power 
            s_deg = r.degree() - d
            s = self.__class__.from_dict( { s_deg : s_coef})
            
            #Update values
            q = s + q
            r = r - s * other
        
        return (q, r)
    
    #TODO: typesetter, talvez definir abstract method na classe PolyOverIntegralDomain
    # e pô-lo lá    
    def gcd(self, other, monic=False):
        """
        Returns the greatest common divisor of self and other with leading 
        coefficient = 1.
        Recall that, if f = g*q + r with deg(r) < deg(g), then
            gcd(f,g) = gcd(g, r),
        and the remainder eventually reaches zero.
        If monic = True then a monic polynomial is returned. Recall that gcd is
        unique up to multiplication by the coefRing's units, which in this case
        are all the non zero elements.
        """
        if self.is_zero() and other.is_zero():
            raise ValueError('gcd(0,0) is not defined.')
        
        f = self.copy()
        g = other.copy()
        
        while not g.is_zero():
            _, rem = f.div_mod(g)
            f = g
            g = rem
        
        result = f
        if monic:
            lead_coef = f.leading_coef()
            result = f._mul_monomial(lead_coef.inverse(), 0)
        
        return result
    
    #TODO: typesetter, talvez definir abstract method na classe PolyOverIntegralDomain
    # e pô-lo lá
    def extended_gcd(self, other, monic=False):
        """
        Extended euclidean algorithm. Returns a tuple (h, s, t) of polynomials 
        satisfying Bezout's identity:
            -> h = gcd(self, other)
            -> h = self * s + other * t
            -> deg(s) < deg(other) - deg(h)
            -> deg(t) < deg(self) - deg(h).
        If monic = True then h is monic and s, t are multiplied by the leading
        coefficient of h.
        """
        r0 = self.copy()
        r1 = other.copy()
        _, r2 = r0.div_mod(r1)
        
        s0 = self.__class__.one()
        s1 = self.__class__.zero()
        
        t0 = self.__class__.zero()
        t1 = self.__class__.one()
        
        while not r2.is_zero():
            q, r2 = r0.div_mod(r1)
            
            s2 = s0 - s1*q
            t2 = t0 - t1*q
            
            r0 = r1
            r1 = r2
            
            s0 = s1
            s1 = s2
            
            t0 = t1
            t1 = t2
            
        if monic:
            lead_coef = r0.leading_coef()
            r0 = r0._mul_monomial(lead_coef.inverse(), 0)
            s0 = s0._mul_monomial(lead_coef.inverse(), 0)
            t0 = t0._mul_monomial(lead_coef.inverse(), 0)
        
        return (r0, s0, t0)
            
            
        
        
        
           


# =============================================================================

class ComRingQuotient( ComRing ):
    
    #Class attribute that says if can reduce to cannonical
    can_reduce = False
    
    @abstractclassmethod
    def is_equivalent(cls, a, b):
        pass
    
    @abstractclassmethod
    def baseRing(cls, other):
        pass
    
    def __init__(self, rep ):
        self.rep = rep
    
    def add(self, other):
        return self.__class__( self.rep + other.rep )
    
    def equals(self, other):
        if self.can_reduce:
            return self.rep == other.rep
        else:
            return self.is_equivalent( self.rep, other.rep )
    
    def mul(self, other):
        return self.__class__( self.rep * other.rep )
    
    @classmethod
    def one(cls):
        if cls.can_reduce:
            return cls( cls.baseRing().one(), reduce = False )
        else:
            return cls( cls.baseRing().one() )
    
    @classmethod
    def zero(cls):
        if cls.can_reduce:
            return cls( cls.baseRing().zero(), reduce = False )
        else:
            return cls( cls.baseRing().zero() )
    
    def symmetric(self):
        return self.__class__( -self.rep )
    
    def is_zero(self):
        if self.can_reduce:
            return self.rep.is_zero()
        else:
            return self.is_equivalent( self.rep, self.zero() )
    
    def is_one(self):
        if self.can_reduce:
            return self.rep.is_one()
        else:
            return self.is_equivalent( self.rep, self.one() )


#TODO: Is it EuclideanDomain?
class ComRingQuotientED( ComRingQuotient ):
    
    """
    Ring that is R/I where R is an EuclideanDomain.
    """
    
    #Overrides
    can_reduce = True
    
    @abstractclassmethod
    def gen(self, other):
        """
        Returns the generator of the Ideal.
        """
        pass
    
    def is_equivalent(self, other):
        return other.isdivisible( self.gen() )
    
    def __init__(self, rep, reduce = True):
        """
        Overrides constructor of ComRingQuotient.
        
        Arguments:
            rep - representative
            reduce - boolean default True - uses rep % gen.
        """
        if reduce:
            self.rep = rep % self.gen()
        else:
            self.rep = rep    

#TODO: Class PolyOverGalois(PolyOverField) with is_irreducible implemented
class PolyOverGalois(PolyOverField):
    
    #TODO:
    def is_irreducible(self):
        pass


#TODO: can then make functions that generate these classes for various gen   
        

