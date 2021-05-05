#pragma once

#include <RcppArmadillo.h>

using namespace arma;

void
kktCheck(uvec& violations,
         uvec& screened,
         const vec& c,
         const uvec& check_set,
         const double lambda)
{
  for (auto&& j : check_set) {
    if (std::abs(c(j)) >= lambda) {
      violations[j] = true;
      screened[j] = true;
    }
  }
}