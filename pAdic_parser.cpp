#include "stdafx.h"
#include <string>
#include <vector>
#include <algorithm>//for sort, see the "order_by_ppower" function
#include "pAdic_parser.h"
#include "pAdic_utilities.h"
#include "pAdic.h"
#include <iostream>

using namespace pAdic_lib;

std::string initial_check_and_clean(const std::string& s, l_int p)
{
	if (s.empty() || s == "0") return "0";
	std::string t{ s };
	for (unsigned int i = 0; i < t.length(); i++)
	{
		//Check for invalid characters
		if (!(t[i] == ' ' || t[i] == 'p' || t[i] == 'P' || t[i] == '+' ||
			t[i] == '-' || t[i] == '*' || t[i] == '^' || t[i] == '(' || t[i] == ')' || isdigit(t[i])))
			throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Syntax_Error);
		//Deleting spaces and brackets
		//To do this deletion, we first check whether (x)y appear, because it cannot turned to be xy. So x)y -> x*y
		if (t[i] == ')' && isdigit(t[i - 1]) && isdigit(t[i + 1])) t.replace(i, 1, "*");
		if (t[i] == ' ' || t[i] == '(' || t[i] == ')') {
			t.replace(i, 1, ""); i--;
		}
	}
	//After deleting brackets we change the occasional +- to -
	for (unsigned int i = 0; i < t.length() - 1; i++)
		if (t[i] == '+' && t[i + 1] == '-') {
			t.replace(i, 1, "");
			t.replace(i, 1, "-");
		}

	//Finally replace all the occurrences of the prime to 'p'
	/*std::string sp{ std::to_string(p) };
	unsigned int sp_len{ sp.length() };
	unsigned int pos{ t.find(sp) };
	unsigned int t_len{ t.length() };
	while (pos != std::string::npos)
	{
		char fdigit{ '.' }; char ldigit{ '.' };
		if (pos != 0) fdigit = t[pos - 1];
		if (pos != t_len - sp_len) ldigit = t[pos + sp_len];
		if (!isdigit(fdigit) && !isdigit(ldigit)) {
			t.erase(pos, sp_len); t.insert(pos, "p");
			t_len = t.length();
		}
		pos = t.find(sp, pos + 1);
	}*/
	return t;
}

std::string initial_check_and_clean_f(const std::string& s)
{
	std::string t{ s };
	for (unsigned int i = 0; i < t.length(); i++)
	{
		//Check for invalid characters
		if (!(t[i] == ' ' || t[i] == '+' ||
			t[i] == '-' || t[i] == '/' || isdigit(t[i])))
			throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Syntax_Error);
		//Deleting spaces and brackets
		//To do this deletion, we first check whether (x)y appear, because it cannot turned to be xy. So x)y -> x*y
		if (t[i] == ' ') { t.replace(i, 1, ""); i--; }
	}
	return t;
}

pAdic parser_f(const std::string& w, l_int p, l_int prec)
{
	l_int n1{ strtol(w.substr(0, w.find('/', 0)).c_str(), nullptr, 10) };
	l_int n2{ strtol(w.substr(w.find('/', 0)+1).c_str(), nullptr, 10) };
	pAdic nn1(n1, p, prec);
	pAdic nn2(n2, p, prec);
	return nn1 / nn2;
}

Vatom parser(const std::string& w, l_int p)
{
	if (w == "")
	{
		return Vatom{ { 0, 0 } };
	}

	//Start parsing: build up a list of consecutive terms (pair of digits and prime powers
	Vatom v{};
	std::string t{ w };
	while (t.length() > 0)
	{
		p_atom pt{};
		signed char sgn{ 1 };
		//the first character is a sign or a digit or 'p'
		if (!(t[0] == '+' || t[0] == '-' || t[0] == 'p' || isdigit(t[0]))) throw pAdic_lib::pAdicError(pAdic_lib::pAdicError::Syntax_Error);
		sgn = (t[0] == '-') ? -1 : sgn = 1;
		if (t[0] == '-' || t[0] == '+') t.replace(0, 1, "");//after stored, sign is deleted
		//if t begins with a digit, it is a coefficient, because the prime is represented by 'p', thanks to initial_check_and_clean
		pt.coeff = sgn;
		if (isdigit(t[0])) {
			int i{};
			while (isdigit(t[i])) i++;
			pt.coeff = sgn * strtol(t.substr(0, i).c_str(), nullptr, 10);
			t.erase(0, i);
			if (t[0] == '*') t.erase(0, 1);
			if (t.length() == 0) {
				pt.ppower = 0;
			}
		}//if (isdigit(t[0]))

		if (t[0] == 'p') {
			if (t[1] == '^') {
				if (t[2] != '-' && t[2] != '+' && !isdigit(t[2])) throw pAdicError(pAdicError::Syntax_Error);
				int i{ (t[2] == '-' || t[2] == '+') ? 3 : 2 };
				while (isdigit(t[i])) i++;
				pt.ppower = strtol(t.substr(2, i - 1).c_str(), nullptr, 10);//includes sign
				t.erase(0, i);
			}
			else {//standalone 'p' without power
				pt.ppower = 1; t.erase(0, 1);
			}
		}//if (t[0] == 'p')
		v.push_back(pt);
	}//while (t.length() > 0)
	normalize_coefficients(v, p);
	order_by_ppower(v);
	return v;
}

void order_by_ppower(Vatom& v)
{
	sort(v.begin(), v.end(), [](p_atom a1, p_atom a2) { return a1.ppower < a2.ppower; });
	return;
}