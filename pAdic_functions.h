#pragma once
#include "pAdic.h"
#include <vector>
#include <utility>

namespace pAdic_lib
{
	pAdic pow(const pAdic&, l_int);
	pAdic exp(const pAdic&);
	pAdic log(const pAdic&);
	pAdic gamma(const pAdic&);
	pAdic exp_AH(const pAdic&);
	//a^b. It is defined only when a is a 1-unit, i.e., when p|(a-1)
	pAdic pow(const pAdic&, const pAdic&);

	using SCF_pairs = std::pair<std::vector<l_int>, std::vector<l_int>>;
	SCF_pairs contfract_Schneider(const pAdic&, const l_int size = -1);
	void print_SCF_pairs(const SCF_pairs& v);
	pAdic SCF_to_pAdic(const SCF_pairs& v, l_int);

	using pAdic_func = pAdic(*) (pAdic);
	pAdic Volkenborn(pAdic_func, l_int, l_int = pAdic::get_global_precision());
};//namespace pAdic_lib
