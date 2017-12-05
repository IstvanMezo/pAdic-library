// pAdic_tryout.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include "pAdic.h"

using namespace std;

//A simple function working with p-adic numbers
pAdic square(pAdic x)
{
	return x * x;
}

int main()
{
	l_int p{ 7 };
	pAdic a(5 * 7 * 7 + 3 * 7 * 7 * 7, 7);
	pAdic b("5*p^0+4*p^1+4*p^2+4*p^3+4*p^4+4*p^5+4*p^6+4*p^7", p);
	cout << a.valuation();

	pAdic::set_global_precision(7);
	pAdic::how_to_print = pAdic_lib::print_style::print_prime;

	cout << a + b << '\n';

	return 0;
}
