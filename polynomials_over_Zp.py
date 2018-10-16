# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 19:31:40 2018

@author: Miguel
"""
from sympy.polys import galoistools as gt
from sympy.ntheory import isprime
from sympy.polys.domains import ZZ


def _remove_leading_zeros(lst):
    """
    Returns a list with deleted last entries of the input that are equal to zero.
    Also, returns an empty list if the input is [0].
    """
    if lst == [] or lst == [0]:
        return []
    else: 
        while lst != [] and lst[-1] == 0:
            lst = lst[:-1]
        return lst
    
def _represents_natural_number(string):
    try:
        i = int(string)
        return i >= 0
    except:
        return False

def make_poly_ring(p):
    """
    Returns a class representing the ring of polynomials over the finite field Z/(p).
    """
    
    assert(isprime(p)), '%d is not a prime number' %p

    class PolynomialsOverZp:
        """
        A polynomial's attribute will be a list with coefficients listed in
        increasing order of degree.
        """
        
        mod = p
        
        def __init__(self, lst):
            # lists the remainders of the integer division by p
            lst = [x % p for x in lst]
            
            # removes any zeros at the end of the list
            lst = _remove_leading_zeros(lst)
            
            self.coefs = lst
        
        @classmethod
        def from_dict(cls, dct):
            """
            In case we want to instantiate a polynomial using a dict.
            Makes it easier to create polynomials such as 1 + x^100.
            Arguments
            We assume the dict has the form {..., power : coefficient, ...}, 
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
        
        def is_zero(self): 
            return self.coefs == []
        
                
        def typesetter(func):
            """
            Decorator.
            Allows binary operations to accept several input types: class instances 
            of polynomials, lists and dicts.
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
                PolynomialsOverZp([1,0,1]).is_equal({'0': 1, '2': 1}.)
            """
            return self.coefs == other.coefs
        
        def __eq__(self, other):
            """
            Overriding the is_equal() method.
            Allows expressions like the following to be evaluated to True:
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
            self_coefs = self.coefs.copy()
            self_coefs.reverse()
            
            other_coefs = other.coefs.copy()
            other_coefs.reverse()
            
            result_coefs = gt.gf_add(self_coefs, other_coefs, p, ZZ)
            result_coefs.reverse()
            result = PolynomialsOverZp(result_coefs)
            return result
        
        def __add__(self, other):
            """
            Overriding. Allows to write expressions like f + g, where f and g
            are class instances.
            """
            return self.add(other)
        
        def neg(self):
            """
            Return the additive inverse in the ring Zp[X].
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
            """
            return self.sub(other)
        
        def to_list(self):
            return self.coefs
        
        def to_dict(self):
            result = {'%d'%i : self.coefs[i] for i in range(
                    len(self.coefs)) if self.coefs[i] != 0}
            return result
        
        def get_degree(self):
            if self.is_zero():
                return '-inf'
            else:
                return len(self.coefs) - 1
        
        #TODO: acertar o output para não aparecerem coisas do tipo
        #                    1*x + 1*x^2.
        def __repr__(self):
            """
            Returns a string with the usual representation of the polynomial.
            """
            if self.is_zero():
                return '0'
            else:
                # turns the coefficients into strings, and 1 into the empty string
                #coefs_to_print = ['' if c == 1 else '%s' %c for c in self.coefs]
                coefs_to_print = self.coefs
                result = ' + '.join(
                        ['%s*x^%s' %(coef, i) for i, coef in enumerate(
                                coefs_to_print) if coef != 0]
                        )
                result = result.replace('x^0', '1')
                return result
    
    
    return PolynomialsOverZp


#==============================================================================
    

Zp_X =  make_poly_ring(2)
h = Zp_X([1, 0, 1, 0, 0, 1, 1, 1])
g = Zp_X([0, 0, 1, 1])

h_coefs = h.coefs
h_coefs.reverse()



g_coefs = g.coefs.copy()
g_coefs.reverse()


teste = gt.gf_add(h_coefs, g_coefs, 3, ZZ)
teste2 = gt.gf_add(g_coefs, h_coefs, 3, ZZ)

dct = {'100' : 1, '0' : 1}
f = Zp_X.from_dict(dct)

f = Zp_X([0,1,2])
g = Zp_X([1,2,0])



methods = ['is_zero', 'to_list', 'to_dict', 'get_degree', '__repr__']

for method in methods:
    result = getattr(h, method)()
    print(method + ':', result)

h.__repr__()

#==============================================================================

print(gt.gf_add([1,1,1], [1,0]))





