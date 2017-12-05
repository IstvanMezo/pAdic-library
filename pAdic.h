#pragma once

#include <string>
#include <vector>
#include <cmath>

namespace pAdic_lib {

	using l_int = long long int;
	using l_double = long double;
	struct p_atom;
	using Vatom = std::vector<p_atom>;

	class pAdicError {
	public:
		enum Error_Types {
			Syntax_Error, NonMatching_primes, IterationError, Division_by_Zero,
			Exp_OutofDomain, Log_OutofDomain, Gamma_OutofDomain, Exp_AH_OutofDomain, Pow_OutofDomain
		};
		pAdicError(Error_Types i) : ErrorType{ i } {}
		Error_Types get_error() { return ErrorType; }
	private:
		Error_Types ErrorType;
	};

	enum class print_style { print_prime, print_p };//with print_prime we print the value of the prime. Otherwise 'p' is printed

	struct p_atom {
		l_int coeff{}, ppower{};
	};

	class pAdic {
	public:
		//constructors
		pAdic(l_int);//initialization by zero, the only parameter is p
		pAdic(const std::string&, l_int, l_int = global_prec);//initialization by a string like "2*7^3-3*7^-5"
		pAdic(const Vatom&, l_int, l_int = global_prec);//initialization by {coefficient, prime power} vector of pairs
		pAdic(l_int, l_int, l_int = global_prec);//initialization by integer
		pAdic(const pAdic& n) = default;//default copy constructor. Note that precision of n is also copied!

		l_int valuation() const { return n[0].ppower; };
		l_double norm() const { return std::pow(p, -valuation()); }

		//operators
		l_int operator[](l_int i) const;//a[i] is the coeff of p^i in the number
		friend pAdic operator+(const pAdic& n1, const pAdic& n2);
		friend pAdic operator+(const l_int n1, const pAdic& n2) { return pAdic(n1,n2.get_p(),n2.get_prec()) + n2; }
		pAdic operator+(const l_int n1) { pAdic ne(n1, p, prec); return *this + ne; }
		pAdic& operator+=(const pAdic& n1) { *this = *this + n1; return *this; }
		pAdic& operator+=(const l_int n1) { *this = *this + n1; return *this; }

		friend pAdic operator-(const pAdic& n1, const pAdic& n2);
		friend pAdic operator-(const l_int n1, const pAdic& n2) { return pAdic(n1, n2.get_p(), n2.get_prec()) - n2; }
		pAdic operator-(const l_int n1) { pAdic ne(n1, p, prec); return *this - ne; }
		pAdic& operator-=(const pAdic& n1) { *this = *this - n1; return *this; }
		pAdic& operator-=(const l_int n1) { *this = *this - n1; return *this; }

		friend pAdic operator*(const pAdic& n1, const pAdic& n2);
		friend pAdic operator*(const l_int n1, const pAdic& n2) { return pAdic(n1, n2.get_p(), n2.get_prec()) * n2; }
		pAdic operator*(const l_int n1) { pAdic sn(n1, p, prec); return *this * sn; }
		pAdic& operator*=(const pAdic& n1) { *this = *this * n1; return *this; }
		pAdic& operator*=(const l_int n1) { *this = *this * n1; return *this; }

		friend pAdic operator/(const pAdic& n1, const pAdic& n2);
		friend pAdic operator/(const l_int n1, const pAdic& n2) { return pAdic(n1, n2.get_p(), n2.get_prec()) / n2; }
		pAdic operator/(const l_int n1) { pAdic sn(n1, p, prec); return *this / sn; }
		pAdic& operator/=(const pAdic& n1) { *this = *this / n1; return *this; }
		pAdic& operator/=(const l_int n1) { *this = *this / n1; return *this; }

		bool operator==(const pAdic& n1);
		bool operator!=(const pAdic& n1) { return !(*this == n1); }

		//auxiliaries
		static print_style how_to_print;
		friend std::ostream& operator<<(std::ostream&, const pAdic&);
		l_int get_prec() const { return prec; }
		l_int get_p() const { return p; }
		void set_prec(l_int = global_prec);
		static void set_global_precision(l_int newprec) { global_prec = newprec; }
		static l_int get_global_precision() { return global_prec; }
		Vatom vectorize() const { return n; }
		bool is_zero() const { if (n.size() == 1 && n[0].coeff == 0) return 1; else return 0; }
		bool is_one() const { if (n.size() == 1 && n[0].coeff == 1 && n[0].ppower == 0) return 1; else return 0; }
		bool is_unit() const { return (valuation() == 0); }
		bool is_1unit() const { return (is_unit() && n[0].coeff == 1); }
	private:
		l_int p;
		Vatom n{};
		l_int prec{ global_prec };
		static l_int global_prec;
	};
}//namespace pAdic
