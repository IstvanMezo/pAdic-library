#include "stdafx.h"
#include "pAdic.h"
#include "pAdic_functions.h"
#include "pAdic_utilities.h"
#include <iostream>
#include <string>

using namespace pAdic_lib;

pAdic pAdic_lib::pow(const pAdic& x, l_int p)
{
	if (p == 0) return pAdic(1, x.get_p(), x.get_prec());
	if (p == 1) return x;
	if (p < 0) return pAdic(1,x.get_p(),x.get_prec()) / pow(x, -p);

	pAdic tmp{ pow(x, p / 2) };
	if (p % 2 == 0) return tmp * tmp;
	else return x * tmp * tmp;
}

pAdic pAdic_lib::exp(const pAdic& x)
{
	if (x.valuation() <= 0) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Exp_OutofDomain);

	pAdic res(1, x.get_p(), x.get_prec());
	pAdic term{ x };
	pAdic pw{ x*x };
	int n{ 2 };
	while( term.valuation() <= x.get_prec() && !term.is_zero() && pw.valuation() <= x.get_prec() && !pw.is_zero())
	{
		res += term;
		term = pw / fact(n);
		n++;
		pw = pow(x, n);
	}
	return res;
}

pAdic pAdic_lib::log(const pAdic& x)
{
	pAdic xm1{ x }; xm1 -= 1;
	if (xm1.valuation() <= 0) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Log_OutofDomain);

	pAdic res(0, x.get_p(), x.get_prec());
	pAdic term{ xm1 };
	pAdic pw{ xm1*xm1 };
	int n{ 2 };
	while (term.valuation() < x.get_prec() && !term.is_zero() && pw.valuation() < x.get_prec() && !pw.is_zero())
	{
		res += term;
		term = pw / (n * (n % 2 == 0 ? -1 : 1));
		n++;
		pw = pow(xm1, n);
	}
	return res;
}

//This is an alternative definition, recursive as the factorial
//l_int int_gamma(l_int x, l_int p)
//{
//	if (x == 1) return -1;
//	else return ( int_gamma(x - 1, p) * ((x - 1) % p == 0 ? -1 : -(x - 1)) );
//}

pAdic pAdic_lib::gamma(const pAdic& x)
{
	if (x.valuation() < 0) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Gamma_OutofDomain);
	if (x.is_zero()) return pAdic(1, x.get_p(), x.get_prec());
	if (x.is_one()) return pAdic(-1, x.get_p(), x.get_prec());

	l_int s{};
	for (auto& e : x.vectorize()) s += e.coeff * int_pow(x.get_p(), e.ppower);
	l_int prod{ 1 };
	for (l_int i = s - 1; i > (s-1) / x.get_p(); i--) prod *= i;
	char sgn{ s % 2 == 0 ? 1 : -1 };
	l_int int_res{ sgn * prod / int_pow(x.get_p(), (s - 1) / x.get_p()) };
	return pAdic(int_res, x.get_p(), x.get_prec());
}

pAdic pAdic_lib::exp_AH(const pAdic& x)
{
	if (x.valuation() <= 0) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Exp_AH_OutofDomain);

	l_int p{ x.get_p() };
	pAdic res{ x };
	pAdic term{ pow(x,p) };	term /= p;
	pAdic pw{ term };
	int i{ 1 };
	while (term.valuation() < x.get_prec() && !term.is_zero() && pw.valuation() < x.get_prec() && !pw.is_zero())
	{
		res += term;
		i++;
		pw = pow(x, int_pow(p, i));
		term = pw / int_pow(p, i);
	}
	return exp(res);

}

pAdic pAdic_lib::pow(const pAdic& a, const pAdic& b)
{
	if (a.get_p() != b.get_p()) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::NonMatching_primes);
	if (!a.is_1unit()) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Pow_OutofDomain);

	if (b.is_one() || a.is_one()) return a;

	return exp(b*log(a));
}

SCF_pairs pAdic_lib::contfract_Schneider(const pAdic& a, const l_int size)
{
	SCF_pairs v{};
	unsigned int k{ size == -1 ? (unsigned int)a.get_prec() : (unsigned int)size };
	v.first.resize(k+1);
	v.second.resize(k+1);
	std::vector<pAdic> al{};
	al.push_back(a);

	for (unsigned int n = 0; n <= k; n++)
	{
		for (unsigned int j = 0; j < a.get_p(); j++)
		{
			//if ((al[n] - j).valuation() > 0)
			if (al[n][0] == j)
			{
				v.first[n] = j;
				break;
			}
		}
		v.second[n] = int_pow(a.get_p(), (al[n] - v.first[n]).valuation());
		//if (al[n] == v.first[n]) break;
		if (n < k) al.push_back(1 / ((al[n] - v.first[n]) / v.second[n]));
	}

	return v;
}

void pAdic_lib::print_SCF_pairs(const SCF_pairs& v)
{
	std::cout << '[';
	std::string comma{ "" };
	for (l_int e : v.first)
	{
		std::cout << comma<< e;
		comma = ",";
	}
	std::cout << "]\n";

	comma = "";
	std::cout << '[';
	for (l_int e : v.second)
	{
		std::cout << comma << e;
		comma = ",";
	}
	std::cout << "]\n";
}

pAdic pAdic_lib::SCF_to_pAdic(const SCF_pairs& v, l_int p)
{
	std::vector<l_int> A(0), B(0);
	A.push_back(v.first[0]);
	B.push_back(1);
	A.push_back(v.first[1] * A[0] + v.second[0] * 1);
	B.push_back(v.first[1] * B[0] + v.second[0] * 0);
	for (unsigned int n = 2; n < v.first.size(); n++)
	{
		A.push_back(v.first[n] * A[n-1] + v.second[n-1] * A[n-2]);
		B.push_back(v.first[n] * B[n-1] + v.second[n-1] * B[n-2]);
	}
	unsigned int m{};
	for (unsigned int n = 1; n < v.first.size(); n++)
	{
		if (B[n] != 0) m = n;
		else break;
	}
	pAdic num( A[m], p ), denum( B[m], p );
	return pAdic(num / denum);
}

pAdic pAdic_lib::Volkenborn(pAdic_func f, l_int p, l_int prec)
{
	pAdic iter{ f(pAdic{0,p,prec}) };
	pAdic old{ iter };//so they are not equal, we enter the while cycle
	
	int n;
	int max{ 3 };
	for (n = 0; n < max; n++)
	{
		old = iter;
		for (l_int j = int_pow(p, n); j < int_pow(p, n + 1); j++)
			iter += f(pAdic(j, p, prec));
		if (old == iter) break;
	}
	return iter / int_pow(p,max);
}