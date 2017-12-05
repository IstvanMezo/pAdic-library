#include "stdafx.h"
#include "pAdic.h"
#include "pAdic_parser.h"
#include "pAdic_utilities.h"
#include <string>
#include <iostream>

using namespace pAdic_lib;

//implementation of the pAdic class

l_int pAdic::global_prec{ 10 };
print_style pAdic::how_to_print{ print_style::print_p };

pAdic::pAdic(l_int pp) : n{ {0,0} }, p{ pp } {};

pAdic::pAdic(const std::string& s, l_int pp, l_int pr) : p{ pp }, prec{ pr }
{
	std::string w{};//Preprocessed string's container
	if (s.find('/') != std::string::npos)
	{
		w = initial_check_and_clean_f(s);
		*this = parser_f(w, p, prec);
		return;
	}

	w = initial_check_and_clean(s, p);
	n = parser(w, p);	
}

pAdic::pAdic(const Vatom& v, l_int pp, l_int pr) : p{ pp }, prec{ pr }
{
	for (auto& e : v) { n.push_back({e.coeff, e.ppower}); }
	normalize_coefficients(n, p);
	order_by_ppower(n);
}

pAdic::pAdic(l_int i, l_int pp, l_int pr) : p{ pp }, prec{ pr }
{
	n = int_to_padic(i, p, prec).n;
}

l_int pAdic::operator[](l_int i) const
{
	for (auto& e : n) { if (e.ppower == i) return e.coeff; }
	return 0;
}

std::ostream& pAdic_lib::operator<<(std::ostream& os, const pAdic& a)
{
	if (a.is_zero()) { os << 0; return os; }
	bool no_first_coeff{ 0 };
	if (a.how_to_print == print_style::print_p)
	{
		for (auto& e : a.n)
		{
			char* sign{ (no_first_coeff) ? "+" : "" };
			if (e.coeff > 0)
			{
				if (e.ppower == 0)
					os << sign << e.coeff;
				else if (e.ppower == 1)
					os << sign << e.coeff << "*p";
				else if (e.ppower >= a.get_prec() - 2)
				{
					if( e.ppower == a.get_prec()-2 ) os << sign << e.coeff << "*p^" << e.ppower;
					os << "+O(p^" << a.get_prec() - 1 << ")";
					return os;
				}
				else os << sign << e.coeff << "*p^" << e.ppower;
			}
			no_first_coeff = 1;
		};
	}
	else
	{
		for (auto& e : a.n)
		{
			char* sign{ (no_first_coeff) ? "+" : "" };
			if (e.coeff > 0)
			{
				if (e.ppower == 0)
					os << sign << e.coeff;
				else if (e.ppower == 1)
					os << sign << e.coeff << "*" << a.p;
				else if (e.ppower >= a.get_prec() - 2)
				{
					if (e.ppower == a.get_prec() - 2) os << sign << e.coeff << "*" << a.p << "^" << e.ppower;
					os << "+O(" << a.p << "^" << a.get_prec() - 1 << ")";
					return os;
				}
				else os << sign << e.coeff << "*" << a.p << "^" << e.ppower;
			}
			no_first_coeff = 1;
		};
	}
	return os;
}

void pAdic_lib::pAdic::set_prec(l_int pr)
{
	for (unsigned int i = 0; i < n.size(); i++)
	{
		if (n[i].ppower > pr) { n.erase(n.begin() + i); i--; }
	}
	prec = pr;
	if (n.size() == 0) n.push_back({ 0,0 });
}

pAdic pAdic_lib::operator+(const pAdic& n1, const pAdic& n2)
{
	//Exception handling
	if (n1.p != n2.p) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::NonMatching_primes);

	l_int min{ l_min( n1.valuation(), n2.valuation() ) };
	l_int max{ l_max( n1.n.back().ppower, n2.n.back().ppower ) };

	Vatom v{};
	bool carry{ 0 };
	l_int i;
	for (i = min; i <= max; i++)
	{
		v.push_back({ (n1[i] + n2[i] + carry) % n1.p,i });
		carry = ((n1[i] + n2[i] + carry) >= n1.p);
	}
	if (carry) { v.push_back({ 1,i }); }

	remove_zeros(v);
	pAdic n(v, n1.p, l_min(n1.get_prec(), n2.get_prec()));
	return n;
}

pAdic pAdic_lib::operator-(const pAdic& n1, const pAdic& n2)
{
	//Exception handling
	if (n1.p != n2.p) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::NonMatching_primes);
	
	pAdic n(n1 + pAdic_minusone(n1.p, l_min(n1.get_prec(), n2.get_prec())) * n2);
	return n;
}

bool pAdic::operator==(const pAdic & n1)
{
	if (n1.p != p) return 0;
	if (n1.n.size() != n.size()) return 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		if (n1.n[i].coeff != n[i].coeff || n1.n[i].ppower != n[i].ppower) return 0;
	}
	return 1;
}

pAdic pAdic_lib::operator*(const pAdic& n1, const pAdic& n2)
{
	//Exception handling
	if (n1.p != n2.p) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::NonMatching_primes);

	Vatom v{};
	for (const auto& e1 : n1.n)
		for (const auto& e2 : n2.n)
		{
			if (e1.coeff * e2.coeff)
			v.push_back({e1.coeff * e2.coeff, e1.ppower + e2.ppower});
		}
	//We do not want to copy the unnecessary vector elements remaining from multiplication, so we delete them
	l_int prec = l_min(n1.get_prec(), n2.get_prec());
	for (unsigned int i = 0; i < v.size(); i++)
	{
		if (v[i].ppower > prec) { v.erase(v.begin() + i); i--; }
	}
	normalize_coefficients(v, n1.p);

	pAdic_lib::pAdic n(v, n1.p, prec);
	//if the result is too big, is must be considered zero up to prec precision.
	if (n.n.size() == 0) n.n.push_back({ 0,0 });
	return n;
}

pAdic pAdic_lib::operator/(const pAdic& n1, const pAdic& n2)
{
	//Exception handling
	if (n2.is_zero()) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Division_by_Zero);
	if (n1.p != n2.p) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::NonMatching_primes);

	l_int alpha = n2.n[0].ppower; 
	l_int prec = l_min(n1.get_prec(), n2.get_prec());
	pAdic n2prime(n2); n2prime.set_prec(prec);
	for (unsigned int i = 0; i < n2prime.n.size(); i++) n2prime.n[i].ppower -= alpha;

	//Newton iteration to find 1/n2
	l_int initpoint{ mod_inverse(n2prime.n[0].coeff, n2.p) };
	pAdic x(0, n2prime.p, prec);//initpoint is simple, between 0..p-1, so no need to call the internal general function int_to_atom, which is, in general, costly
	x.n[0].coeff = initpoint;
	x.n[0].ppower = 0;

	pAdic two(0, n2prime.p, prec);
	if (n2prime.p == 2) two.n[0] = { 1,1 }; else two.n[0] = { 2,0 };
	pAdic minusone(pAdic_minusone(n1.p, prec));
	pAdic minusone_x_n2prime = minusone * n2prime;

	pAdic oldx{ 0,n1.p,prec };
	int iternum{ 1 };
	const char maxiternum = 20;
	do
	{
		oldx = x;
		x = (two + minusone_x_n2prime * x) * x;
		iternum++;
	} while (iternum < maxiternum && x != oldx );
	if (iternum == maxiternum) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::IterationError);
	
	pAdic res(n1 * x);
	for (auto& e : res.n) e.ppower -= alpha;
	//if the divisor is too big, the result is zero up to prec precision. To assure that the vector n is not empty, we add {0,0}
	if (res.n.size() == 0) res.n.push_back({ 0,0 });
	return res;
}