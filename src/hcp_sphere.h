#ifndef HCP_H
#define HCP_H

#include "utils.h"
#include "sampling.h"
#include "geometric_algorithms.h"
#include "sort.h"
#include "hilbert_curve.h"

#include <iostream>
#include <set>

using iterator = std::vector<ptrdiff_t>::iterator;


ints curve_order(const vec2s& X);
vecs general_plan(const double* a_weight, int n, const double* b_weight, int m);
vecs General_Plan(const scalars& X, const scalars& Y);




#endif // HCP_H
