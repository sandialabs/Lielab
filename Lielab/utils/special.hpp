#ifndef LIELAB_UTILS_SPECIAL_HPP
#define LIELAB_UTILS_SPECIAL_HPP

#include <vector>

namespace Lielab::utils
{

double bernoulli(const int num, const int sign = -1);
int sign(const double x);
double sinc(const double x);

}

#endif
