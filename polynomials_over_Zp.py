# -*- coding: utf-8 -*-
"""
"""

# =============================================================================

from sympy.polys import galoistools as gt
from sympy.ntheory import isprime
from sympy.polys.domains import ZZ

# =============================================================================

#teste

def _remove_trailing_zeros(lst):
    """
    Removes any zeros at the end of the list.
    """
    k=0
    for k, value in enumerate( lst[::-1] ):
        if value != 0:
            break
    lst_no_trailing_zeroes = lst if k == 0 else lst[:-k]
    return lst_no_trailing_zeroes

def _represents_natural_number(string):
    """
    Returns True if the input string is an integer >= 0, False otherwise.
    """
    try:
        i = int(string)
        return i >= 0
    except:
        return False

def make_poly_ring(p):
    """
    Arguments:
        p - prime number >= 2.
    Returns a class representing the ring of polynomials over the finite 
    field Z/(p) =: Zp.
    """
    
    assert(isprime(p)), '%d is not a prime number' %p

    class PolynomialsOverZp:
        """
        A polynomial's attribute will be a list with coefficients listed in
        increasing order of degree.
        """
        
        #class attribute
        mod_p = p
        
        def __init__(self, lst):
                        
            # lists the remainders of the integer division by p
            lst = [x % p for x in lst]
            
            # removes any zeros at the end of the list
            lst = _remove_trailing_zeros(lst)
            
            self.coefs = lst
        
        @classmethod
        def from_dict(cls, dct):
            """
            In case we want to instantiate a polynomial using a dict.
            Makes it easier to create polynomials such as 1 + x^100.
            It is assumed the dict has the form {..., power : coefficient, ...}, 
            where type(power) = string and type(coefficient) = int.
            Coefficients are assumed to be zero if they have no corresponding key.
            """
            if len(dct) == 0:
                result = cls([])
            else:
                assert(
                        all([_represents_natural_number(key) for key in dct.keys()])
                        ), 'Invalid keys in dict'
                
                powers = [int(x) for x in dct.keys()]
                degree = max(powers)
                
                # An m-degree polynomial has m+1 coefficients
                lst = [0]*(degree+1)
                
                # fills in the coefficients
                for i in powers:
                    lst[i] = dct[str(i)]
                result = cls(lst)
            
            return result
        
        @classmethod
        def from_int(cls, int_):
            """
            Allows the instatiation of a constant polynomial from an int type,
            as in
                f = PolynomialsOverZp(2).
            Creates an element of the ring Zp viewed as a subring of Zp[X].
            """
            assert(type(int_) == int), 'Input is not of type int.'
            return cls([int_])
        
        def is_zero(self): 
            return self.coefs == []
                
        def typesetter(func):
            """
            Decorator.
            Allows binary operations to accept several input types: class instances 
            of polynomials, lists, dicts and ints.
            """
            def new_func(self, other):
                
                # in this case there is nothing to do, the default methods assume
                # that the input is a class instance
                if isinstance(other, PolynomialsOverZp):
                    result = func(self, other)
                # instantiates a polynomial from the list
                elif type(other) == list:
                    other = PolynomialsOverZp(other)
                    result = func(self, other)
                # instantiates a polynomial from the dict
                elif type(other) == dict:
                    other = PolynomialsOverZp.from_dict(other)
                    result = func(self, other)
                #TODO: ver se não há problemas com o tipo int nos vários métodos.
                #audaxH - não percebi
                # instatiates a constant polynomial
                elif type(other) == int:
                    other = PolynomialsOverZp([other])
                    result = func(self, other)
                else:
                    error_msg = 'Not a valid polynomial.'
                    raise TypeError(error_msg)
                return result
            
            return new_func
        
        @typesetter
        def is_equal(self, other):
            """
            Returns True if self = other in Zp[X], False otherwise. The equality
            is coefficient-wise.
            The method assumes that the argument "other" is a class instance, 
            but the typesetter decorator allows the following expressions to be
            evaluated to True:
                PolynomialsOverZp([1, 0, 1]).is_equal([1, 0, 1])
            and
                PolynomialsOverZp([1,0,1]).is_equal({'0': 1, '2': 1}).
            """
            return self.coefs == other.coefs
        
        def __eq__(self, other):
            """
            Overriding the is_equal() method.
            Allows expressions like the following to be evaluated to True:
                Zp_X([1, 0, 1]) == Zp_X.from_dict({'0':1, '2':1}).
            The decorator in the is_equal() method also allows 
                PolynomialsOverZp([1, 0, 1]) == [1, 0, 1]
            and
                PolynomialsOverZp([1,0,1]) == {'0': 1, '2': 1}.
            """
            return self.is_equal(other)
        
        @typesetter
        def add(self, other):
            """
            Polynomial sum in Zp[x].
            The method assumes that the argument "other" is a class instance, 
            but the typesetter decorator allows the following:
                PolynomialsOverZp([1, 0, 1]).add([0, 1, 0, 1])
            and
                PolynomialsOverZp([1, 0, 1]).add({'1' : 1, '3' : 1}).
            """
            # the copy assures that the self instance is not changed.
            self_coefs = self.coefs.copy()
            
            # the functions that were implemented in sympy.polys.galoistools
            # assume that the coefficients are in decreasing order of degree,
            # unlike this class.
            self_coefs.reverse()
            
            # same as above
            other_coefs = other.coefs.copy()
            other_coefs.reverse()
            
            result_coefs = gt.gf_add(self_coefs, other_coefs, p, ZZ)
            
            # lists the coefficients in increasing order of degree
            result_coefs.reverse()
            
            result = PolynomialsOverZp(result_coefs)
            return result
        
        def __add__(self, other):
            """
            Overriding. Allows to write expressions like f + g, where f and g
            are class instances.
            Also allows things like
                f + [1, 0, 1]
            and
                f + {'100' : 1},
            because of the use of the typesetter decorator in the add() method.
            """
            return self.add(other)
        
        def neg(self):
            """
            Return the additive inverse of "self" in the ring Zp[X].
            Remark: not necessary if p=2, since h = -h in Z2[X].
            """
            return PolynomialsOverZp([(-n) % p for n in self.coefs])
        
        def __neg__(self):
            """
            Overriding. Allows one to write -h for some class instance h.
            """
            return self.neg()
        
        @typesetter
        def sub(self, other):
            """
            Polynomial difference in Zp[X].
            The method assumes that the argument "other" is a class instance, 
            but the typesetter decorator allows the following:
                PolynomialsOverZp([1, 0, 1]).sub([0, 1, 0, 1])
            and
                PolynomialsOverZp([1, 0, 1]).sub({'1' : 1, '3' : 1}).
            """
            return self.add(other.neg())
        
        def __sub__(self, other):
            """
            Overriding. Allows to write expressions like f - g, where f and g
            are class instances.
            Also allows things like
                f - [1, 0, 1]
            and
                f - {'100' : 1},
            because of the use of the typesetter decorator in the sub() method.
            """
            return self.sub(other)
        
        #TODO: docstring
        @typesetter
        def mul(self, other):
            """
            Polynomial multiplication in Zp[X].
            """
            # the copy assures that the self instance is not changed
            self_coefs = self.coefs.copy()
            # the functions that were implemented in sympy.polys.galoistools
            # assume that the coefficients are in decreasing order of degree,
            # unlike this class
            self_coefs.reverse()
            
            # same as above
            other_coefs = other.coefs.copy()
            other_coefs.reverse()
            
            result_coefs = gt.gf_mul(self_coefs, other_coefs, p, ZZ)
            
            # lists the coefficients in increasing order of degree
            result_coefs.reverse()
            
            result = PolynomialsOverZp(result_coefs)
            return result
        
        def __mul__(self, other):
            """
            Overriding. Allows to write expressions like f*g, where f and g
            are class instances.
            Also allows things like
                f * [1, 0, 1]
            and
                f * {'100' : 1},
            because of the use of the typesetter decorator in the mul() method.
            """
            return self.mul(other)
       
        #TODO: testar
        #TODO: doctring
        @typesetter
        def div_mod(self, other):
            """
            Division with remainder in Zp[X].
            Recall that, since Zp is a field, then Zp[X] is an euclidean domain
            and thus this division is possible whenever "other" != 0.
            """
            # the copy assures that the self instance is not changed
            self_coefs = self.coefs.copy()
            # the functions that were implemented in sympy.polys.galoistools
            # assume that the coefficients are in decreasing order of degree,
            # unlike this class
            self_coefs.reverse()
            
            # same as above
            other_coefs = other.coefs.copy()
            other_coefs.reverse()
            
            quot_coefs, remainder_coefs = gt.gf_div(self_coefs, other_coefs, 
                                                    p, ZZ)
            
            # lists the coefficients in increasing order of degree
            quot_coefs.reverse()
            remainder_coefs.reverse()
            
            quot = PolynomialsOverZp(quot_coefs)
            remainder = PolynomialsOverZp(remainder_coefs)
            
            return quot, remainder
        
        #TODO: docstring
        def quotient(self, other):
            """
            Returns the quotient of the euclidean division in Zp[X].
            """
            result, _ = self.div_mod(other)
            return result
        
        #TODO: docstring
        def __floordiv__(self, other):
            """
            Overriding. Allows to write f // g.
            """
            result = self.quotient(other)
            return result
        
        #TODO: testar
        #TODO: docstring
        def mod(self, other):
            """
            Returns the remainder of the euclidean division in Zp[X].
            """
            _, result = self.div_mod(other)
            return result
        
        def __mod__(self, other):
            """
            Overriding. Allows to write f % g.
            """
            return self.mod(other)
        
        #TODO: testar
        #TODO: docstring
        def is_irreducible(self):
            """
            Tests irreducibility.
            """
            if self.is_zero():
                return False
            else:
                # the copy assures that the self instance is not changed
                self_coefs = self.coefs.copy()
                # the functions that were implemented in sympy.polys.galoistools
                # assume that the coefficients are in decreasing order of degree,
                # unlike this class
                self_coefs.reverse()
                
                result = gt.gf_irreducible_p(self_coefs, p, ZZ)
                return result
            
        def to_list(self):
            return self.coefs
        
        def to_dict(self):
            result = {'%d'%i : self.coefs[i] for i in range(
                    len(self.coefs)
                    ) if self.coefs[i] != 0}
            return result
        
        def degree(self):
            if self.is_zero():
                return -1
            else:
                return len(self.coefs) - 1
        
        def __repr__(self):
            """
            Returns a string with the usual representation of the polynomial.
            It is in LaTeX compatible form.
            """
            
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
                    if self.coefs[0] != 0:
                        poly_str_lst.append( str(self.coefs[0]) )
                    first_coef = self.coefs[1]
                    if first_coef != 0:
                        first_coef_ = '' if first_coef == 1 else first_coef
                        poly_str_lst.append( '%sx' % first_coef_ )
                    if deg == 1:
                        poly_str = ' + '.join( poly_str_lst )
                    else:
                        #adds powers to the rest, and hides coef if it is 1
                        for power, coef in zip( range(2, deg+1), self.coefs[2:] ):
                            if coef != 0:
                                coef_ = '' if coef == 1 else coef
                                poly_str_lst.append( '%sx^{%s}' % ( coef_, power) )
                        poly_str = ' + '.join( poly_str_lst )
                    return poly_str
    
    
    return PolynomialsOverZp


#==============================================================================


