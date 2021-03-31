#pragma once

#include <RcppArmadillo.h>

using namespace arma;

inline double
squaredNorm(const vec& x)
{
  return std::pow(norm(x), 2);
}

template<typename T>
inline int
signum(T val)
{
  return (T(0) < val) - (val < T(0));
}

template<typename T, typename S>
inline bool
contains(const T& x, const S& what)
{
  return std::find(x.begin(), x.end(), what) != x.end();
}

inline uvec
setUnion(const uvec& a, const uvec& b)
{
  std::vector<unsigned> out;
  out.reserve(a.n_elem + b.n_elem);

  std::set_union(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  out.shrink_to_fit();

  return conv_to<uvec>::from(out);
}

inline uvec
setDiff(const uvec& a, const uvec& b)
{
  std::vector<uword> out;
  out.reserve(a.n_elem);

  std::set_difference(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  out.shrink_to_fit();

  return conv_to<uvec>::from(out);
}

inline uvec
setIntersect(const uvec& a, const uvec& b)
{
  std::vector<uword> out;
  out.reserve(std::min(a.n_elem, b.n_elem));

  std::set_intersection(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  out.shrink_to_fit();

  return conv_to<uvec>::from(out);
}

// set intersection that retains permutation in `a`
inline uvec
safeSetIntersect(const uvec& a, const uvec& b)
{
  std::vector<uword> out;
  out.reserve(std::min(a.n_elem, b.n_elem));

  for (auto&& a_i : a) {
    if (contains(b, a_i)) {
      out.emplace_back(a_i);
    }
  }

  out.shrink_to_fit();

  return conv_to<uvec>::from(out);
}

// set difference that retains permutation in `a`
inline uvec
safeSetDiff(const uvec& a, const uvec& b)
{
  std::vector<uword> out;
  out.reserve(a.n_elem);

  for (auto&& a_i : a) {
    if (!contains(b, a_i)) {
      out.emplace_back(a_i);
    }
  }

  out.shrink_to_fit();

  return conv_to<uvec>::from(out);
}

inline vec
matTransposeMultiply(const mat& A,
                     const vec& b,
                     const vec& offset,
                     const bool standardize)
{
  return A.t() * b;
}

inline vec
matTransposeMultiply(const sp_mat& A,
                     const vec& b,
                     const vec& offset,
                     const bool standardize)
{
  vec Atb = A.t() * b;

  if (standardize) {
    Atb -= offset * sum(b);
  }

  return Atb;
}

inline vec
matTransposeMultiply(const mat& A,
                     const vec& b,
                     const uvec& ind,
                     const vec& offset,
                     const bool standardize)
{
  return A.cols(ind).t() * b;
}

inline vec
matTransposeMultiply(const sp_mat& A,
                     const vec& b,
                     const uvec& ind,
                     const vec& offset,
                     const bool standardize)
{
  vec Atb = A.cols(ind).t() * b;

  if (standardize) {
    Atb -= offset * sum(b);
  }

  return Atb;
}

double
getSparsity(const mat& X)
{
  return 1;
}

double
getSparsity(const sp_mat& X)
{
  return 1 - X.n_nonzero / X.n_elem;
}
