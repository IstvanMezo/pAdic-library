#include "stdafx.h"
#include "pAdic.h"
#include "pAdic_utilities.h"
#include <vector>
#include <math.h>
#include <limits.h>

using namespace pAdic_lib;

pAdic pAdic_minusone(l_int p, l_int prec)
{
	Vatom v{};
	for (int i = 0; i <= prec; i++) v.push_back({ p - 1,i });
	pAdic n(v, p, prec);
	return n;
}

l_int int_pow(l_int x, l_int p)
{
	if (p == 0) return 1;
	if (p == 1) return x;

	l_int tmp = int_pow(x, p / 2);
	if (p % 2 == 0) return tmp * tmp;
	else return x * tmp * tmp;
}

l_int mod_inverse(l_int a, l_int p)
{
	if (a == 1) return 1;
	l_int t, nt, r, nr, q, tmp;
	if (a < 0) a = p - (-a % p);
	t = 0;  nt = 1;  r = p;  nr = a % p;
	while (nr != 0) {
		q = r / nr;
		tmp = nt;  nt = t - q*nt;  t = tmp;
		tmp = nr;  nr = r - q*nr;  r = tmp;
	}
	if (r > 1) return -1;  /* No inverse */
	if (t < 0) t += p;
	return t;
}

pAdic int_to_padic(l_int n, l_int p, l_int prec)
{
	pAdic a{ int_to_atom(n, p, prec), p, prec };
	return a;
}

Vatom int_to_atom(l_int n, l_int p, l_int prec)
{
	if (n == 0) return Vatom{ { 0, 0 } };
	signed char sgn{ (n < 0) ? -1 : 1 };
	if (sgn == -1) n = -n;

	Vatom v{};
	n %= int_pow(p, prec);//cut off the part greater than the default precision with which the return value will be initialized

	for (auto i = prec; i >= 0; i--)
	{
		auto pp{ int_pow(p, i) };
		auto ac{ n / pp };
		if (ac) {
			v.push_back({ ac,i });
			n -= ac * int_pow(p, i);
		}
	}

	if (sgn == -1)
	{
		Vatom all_pminus1{};
		for (l_int i = 0; i <= prec; i++) all_pminus1.push_back({ p - 1,i });
		for (auto& e : v)
		{
			all_pminus1.erase(all_pminus1.begin() + (int)e.ppower);//(int)e.ppower is to avoid C4244 warning
			all_pminus1.push_back({ p - 1 - e.coeff, e.ppower });
		}
		all_pminus1.push_back({ 1,0 });//the above line is the additive inverse of -n-1, so we still must add 1 = 1*p^0
		v = all_pminus1;
	}
	return v;
}

void normalize_coefficients(pAdic_lib::Vatom& v, l_int p)
{
	Vatom aux{};
	unsigned int i{ 0 };
	//substitute the p-adic equivalents of the coefficients out of [0,p-1]
	while (i < v.size())
	{
		if (v[i].coeff < 0 || v[i].coeff >= p)
		{
			Vatom tmp{ int_to_atom(v[i].coeff, p) };
			for (auto& e : tmp) e.ppower += v[i].ppower;
			aux.insert(aux.end(), tmp.begin(), tmp.end());
			v.erase(v.begin() + i);
		}
		else i++;
	}
	v.insert(v.end(), aux.begin(), aux.end());
	//sum the same p powers by coefficients if there are any
	bool was_match{ 0 };
	for (unsigned int i = 0; i < v.size(); i++)
	{
		for (unsigned int j = i + 1; j < v.size(); j++)
		{
			if (v[j].ppower == v[i].ppower)
			{
				was_match = 1;
				v[i].coeff += v[j].coeff;
				v.erase(v.begin() + j);
				j--;
			}
		}
	}
	if (was_match) normalize_coefficients(v, p);
}

void remove_zeros(pAdic_lib::Vatom& v)
{
	for (unsigned int i = 0; i < v.size(); i++)
		if (v[i].coeff == 0) v.erase(v.begin() + i);
}