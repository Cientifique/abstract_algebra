#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

# =============================================================================

from concrete_algebraic import Integer, PolyOverZ, PolyOverQ, Rational, \
    makeZpX, makeZp

from polynomials_over_Zp import make_poly_ring

#generates random polynomials for testing
import random
import time
# =============================================================================

PRINT_POLY_OVER=False
if PRINT_POLY_OVER:
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


# =============================================================================


#GENERATE RANDOM POLYS

PRIME_P = 11
MAX_DEG = 30
NPOLY = 2
NTESTS = 100

random_coefs = []
for k in range( NPOLY ):
    deg = random.randint(0, MAX_DEG)
    #non-zero poly
    leading_coef = [ random.randint( 1, PRIME_P-1 ) ]
    coefs = [ random.randint( 0, PRIME_P-1 ) for x in range(deg) ]
    random_coefs.append( coefs + leading_coef )
    
# =============================================================================
       
def make_random_poly_list( random_coefs, PolyMaker, from_int_coefs = False ):
    random_poly_lst = []
    for coefs_ in random_coefs:
        if from_int_coefs:
            random_poly_lst.append( PolyMaker.from_int_coefs( coefs_ ) )
        else:
            random_poly_lst.append( PolyMaker( coefs_) )
    return random_poly_lst

# =============================================================================

#using galois tools
ZpX_Gal = make_poly_ring(PRIME_P)
#not using galois tools
ZpX = makeZpX( PRIME_P )

random_pols_Gal = make_random_poly_list(random_coefs, ZpX_Gal )
random_pols = make_random_poly_list(random_coefs, ZpX, True)

# =============================================================================

#to debug
random_poly_lst = random_pols

def poly_tester( random_poly_lst ):
    
    times_dict = {'euclidean_div' : [], 'com_product' : [],
                  'com_addition' : []}
    
    #test with pairs of these
    for k in range(NTESTS):
        
        g = random.choice( random_poly_lst )
        h = random.choice( random_poly_lst )
        
        error_msg = (
                '\n' + '-' * 50 + 'g =' + str( g ) + 'h =' + str( h) + '-' * 50
        )
        
        
        tic = time.process_time()
        if not (h == (h//g)*g + (h%g) ):
            error_ = '\n' + 'failed h == (h//g)*g + (h%g)' + error_msg
            raise ValueError( error_ )
        times_dict['euclidean_div'].append( time.process_time()-tic )
        
        tic = time.process_time()
        if not ( h+g == g+h ):
            error_ = '\n' + 'failed h+g == g+h' + error_msg
            raise ValueError( error_ )
        times_dict['com_addition'].append( time.process_time()-tic )
        
        tic = time.process_time()
        if not (h*g == g*h):
            error_ = 'failed h*g == g*h' + error_msg
            raise ValueError( error_ ) 
        times_dict['com_product'].append( time.process_time()-tic )
            
    return times_dict
#        print('-' * 100)
#        msg_ = 'Test irreducible'
#        print(msg_ + ( 96 -len(msg_) ) * ' ', 
#              g.is_irreducible(), h.is_irreducible()
#              )
#        print('\n')

print('\n\nGALOIS\n\n')
times_Gal = poly_tester(random_pols_Gal)
print('\n\nPURE\n\n')
times_pure = poly_tester(random_pols)

for key in times_pure.keys():
    gal_ = times_Gal[key]
    pure_ = times_pure[key]
    
    print( 'Times for ', key, '\n\n' )
    for name, lst in zip( ('gal', 'pure'), [gal_, pure_] ):
        print( name)
        sum_ = sum(lst) * 1000
        n_ = len(lst)
        print(' '*4, 'n    -> %d tests' % n_ )
        print(' '*4, 'sum  -> %.2f ms' % sum_ )
        print(' '*4, 'mean -> %.2f ms' % (sum_/n_ ) )
        print('\n')

# =============================================================================

#Compare

METHODS = ['__add__', '__mul__', '__mod__', '__floordiv__' ]
    
def comparer( g, h, g_Gal, h_Gal):
    for method in METHODS:
        result = getattr(g, method)( h )
        result_gal = getattr(g_Gal, method)( h_Gal )
        error_msg = '\n'.join( [ str(g), str(h), str(g_Gal), str(h_Gal), method])
        assert( str(result_gal) == str(result) ), error_msg
    
    
for x in range( NTESTS ):
    idx_g = random.randint(0, len(random_pols) - 1) 
    idx_h = random.randint(0, len(random_pols)  - 1) 
    
    g = random_pols[ idx_g ]
    h = random_pols[ idx_h ]
    
    g_Gal = random_pols_Gal[ idx_g ]
    h_Gal = random_pols_Gal[ idx_h ]
    
    comparer( g, h, g_Gal, h_Gal )

# =============================================================================

#Correu sem erros!
print('Há aí o Galois Tools que é o CESSO, mas o CESSO agora é outro')


print('É mais lento mas agora é uma questão de otimizar...',
      'o produto de polinómios pode ser otimizado')

# =============================================================================


#
#print('\n\n', '-' * 100, '\n\n')
#h = Zp_X.from_dict({'31' : 2, '4' : 2, '11' : 1, '0' : 1})
#g = Zp_X([1, 0, 2, 0])
#
#print('h =', h)
#print('g =', g)
#
#quot, rem = h.div_mod(g)
#
#h_coefs = h.coefs
#h_coefs.reverse()
#   
#
#g_coefs = g.coefs.copy()
#g_coefs.reverse()
#
#
#teste = gt.gf_add(h_coefs, g_coefs, 3, ZZ)
#teste2 = gt.gf_add(g_coefs, h_coefs, 3, ZZ)
#
#dct = {'100' : 1, '0' : 1}
#f = Zp_X.from_dict(dct)
#
#f = Zp_X([0,1,2])
#g = Zp_X([1,2,0])
#
#
#
#methods = ['is_zero', 'to_list', 'to_dict', 'degree', '__repr__']
#
#for method in methods:
#    result = getattr(h, method)()
#    print(method + ':', result)
