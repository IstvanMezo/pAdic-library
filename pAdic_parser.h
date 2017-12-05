#pragma once

#include <string>
#include <vector>
#include "pAdic.h"

using namespace pAdic_lib;

//checks whether s contains non-permitted characters, and cleans it from spaces
//then replaces all the occasions of the number p to 'p'
std::string initial_check_and_clean(const std::string&, l_int);

//checks whether s contains non-permitted characters, and cleans it from spaces
//whe the initializer string is of the form "a/b"
std::string initial_check_and_clean_f(const std::string&);


Vatom parser(const std::string&, l_int);

pAdic parser_f(const std::string&, l_int, l_int);

//parser leaves vector<atom> unordered
//order_by_ppower orders the vector increasing by prime powers
void order_by_ppower(Vatom&);