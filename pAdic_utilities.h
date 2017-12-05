#pragma once

#include "pAdic.h"

using namespace pAdic_lib;

//minus one. This explicit function is quicker than the general int_to_padic
pAdic pAdic_minusone(l_int, l_int);//prime and precision must be specified

//quick integer power x^p function
l_int int_pow(l_int x, l_int p);

//factorial
inline constexpr l_int fact(l_int n)
{
	if (n == 1) return 1; else return n*fact(n - 1);
}

//Solves ax==1 (mod p)
//Code taken from https://rosettacode.org/wiki/Modular_inverse#C
l_int mod_inverse(l_int, l_int);

inline l_int l_min(l_int x, l_int y) { return (x < y) ? x : y; }
inline l_int l_max(l_int x, l_int y) { return (x > y) ? x : y; }

//develops n into a padic number
pAdic int_to_padic(l_int n, l_int p, l_int prec = pAdic::get_global_precision());

//develops n into padic form and gives back the p_atom vector
Vatom int_to_atom(l_int n, l_int p, l_int prec = pAdic::get_global_precision());

//recalculates the vector so that all the .coeffs are in [0,p-1]
void normalize_coefficients(pAdic_lib::Vatom& v, l_int p);

//removes those vector elements for which the coeff is zero
//in some cases (like (1 + p^2) - 1) Vatom elements become zero.
//Thus unneccessary elements does not slow computations, and valuation() works simpler
void remove_zeros(pAdic_lib::Vatom& v);