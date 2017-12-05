// pAdic_tryout.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include "pAdic.h"
#include "pAdic_parser.h"
#include "pAdic_utilities.h"
#include "pAdic_functions.h"

using namespace std;

pAdic square(pAdic x)
{
	return x * x;
}


pAdic expon(pAdic x)
{
	static pAdic base{ pAdic("1+p",3) };//the function is called very often from Volkenborn, static base is better
	if (x.is_zero()) return pAdic(1, 3);
	return pow(base,x);
}

int main()
{
	l_int p{ 17 };
	pAdic a(5 * 7 * 7 + 3 * 7 * 7 * 7, 7);
	cout << a.valuation();

	//pAdic::set_global_precision(7);
	//pAdic::how_to_print = pAdic_lib::print_style::print_prime;

	//cout << Volkenborn(square,p) << '\n';
	//cout << pAdic("1/6", p) << '\n';
	//cout << Volkenborn(expon,p) << '\n';
	//cout << pAdic("p^-1", p) * log(pAdic("1+p", p)) << '\n';


	//pAdic c("5 * p ^ 0 + 4 * p ^ 1 + 4 * p ^ 2 + 4 * p ^ 3 + 4 * p ^ 4 + 4 * p ^ 5 + 4 * p ^ 6 + 7 * p ^ 9", p);
	//pAdic c("5 * p ^ 0 + 4 * p ^ 1 + 4 * p ^ 2 + 4 * p ^ 3 + 4 * p ^ 4 + 4 * p ^ 5 + 4 * p ^ 6 + 4 * p ^ 7", p);

	//p = 17;
	//pAdic_lib::pAdic c("-318p^(5)+14*17 ^ (5) + (24) *      17^3+ 28*17^(2)+2+58*p^2+97*17^3", p);
	//cout << c << '\n';

	//pAdic_lib::pAdic d({ {-318,5},{14, 5}, {24, 3 }, { 28,2 }, { 2,0 }, { 58, 2 }, { 97, 3 }}, p);
	//cout << d << '\n';

	//int p{ 17 };
	//pAdic_lib::pAdic b("+(-4117)17^(31)-318p^(-1956)+-3*17 ^ (-5) + (4) *      (17)^-3+ (-2)17^(-2)+2+17", p);
	//b.how_to_print = pAdic_lib::print_style::print_p;
	//cout << b << '\n';
	//cout << b.valuation() << '\n';

	return 0;
}