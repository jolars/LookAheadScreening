#pragma once

#include "utils.h"

double
prox(const double x, const double lambda)
{
  return signum(x) * std::max(std::abs(x) - lambda, 0.0);
}
