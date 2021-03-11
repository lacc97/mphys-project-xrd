#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstdint>

#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "allocator.hpp"
#include "complex.hpp"

#ifndef WITH_LOW_PRECISION
//#  ifdef VERY_HIGH_PRECISION
//#    include <boost/multiprecision/gmp.hpp>
//
////very high precision
//typedef boost::multiprecision::mpz_int sint_t;
//typedef boost::multiprecision::mpz_int uint_t;
//typedef boost::multiprecision::mpf_float_100 real_t;
//#  else
// high precision   (for actual results)
typedef __int128 sint_t;
typedef unsigned __int128 uint_t;
typedef long double real_t;
//#  endif
#else
// low precision   (for quick prototyping and testing)
typedef int64_t sint_t;
typedef uint64_t uint_t;
typedef double real_t;
#endif

typedef math::complex cplx_t;
constexpr auto k_i = cplx_t{0, 1};

constexpr auto n_dynamic = Eigen::Dynamic;

/// Dynamic length multidimensional arrays
template <typename Numeric>
using mdata_t = Eigen::Array<Numeric, n_dynamic, n_dynamic>;
typedef mdata_t<real_t> rmdata_t;
typedef mdata_t<cplx_t> cmdata_t;

/// Dynamic length arrays
template <typename Numeric>
using data_t = Eigen::Array<Numeric, n_dynamic, 1>;
template <typename Numeric>
using data_span_t = Eigen::Map<data_t<Numeric>>;
template <typename Numeric>
using data_view_t = Eigen::Map<const data_t<Numeric>>;
using idata_t = data_t<sint_t>;
using idata_span_t = data_span_t<sint_t>;
using idata_view_t = data_view_t<sint_t>;
using rdata_t = data_t<real_t>;
using rdata_span_t = data_span_t<real_t>;
using rdata_view_t = data_view_t<real_t>;
using cdata_t = data_t<cplx_t>;
using cdata_span_t = data_span_t<cplx_t>;
using cdata_view_t = data_view_t<cplx_t>;

//  hacky?
namespace Eigen {
  template <typename Numeric>
  inline auto begin(data_t<Numeric>& d) noexcept {
    return d.data();
  }
  template <typename Numeric>
  inline auto begin(const data_t<Numeric>& d) noexcept {
    return d.data();
  }
  template <typename Numeric>
  inline auto begin(data_span_t<Numeric>& d) noexcept {
    return d.data();
  }
  template <typename Numeric>
  inline auto begin(const data_span_t<Numeric>& d) noexcept {
    return d.data();
  }
  template <typename Numeric>
  inline auto begin(data_view_t<Numeric>& d) noexcept {
    return d.data();
  }
  template <typename Numeric>
  inline auto begin(const data_view_t<Numeric>& d) noexcept {
    return d.data();
  }
  template <typename Numeric>
  inline auto end(data_t<Numeric>& d) noexcept {
    return d.data() + d.size();
  }
  template <typename Numeric>
  inline auto end(const data_t<Numeric>& d) noexcept {
    return d.data() + d.size();
  }
  template <typename Numeric>
  inline auto end(data_span_t<Numeric>& d) noexcept {
    return d.data() + d.size();
  }
  template <typename Numeric>
  inline auto end(const data_span_t<Numeric>& d) noexcept {
    return d.data() + d.size();
  }
  template <typename Numeric>
  inline auto end(data_view_t<Numeric>& d) noexcept {
    return d.data() + d.size();
  }
  template <typename Numeric>
  inline auto end(const data_view_t<Numeric>& d) noexcept {
    return d.data() + d.size();
  }
}    // namespace Eigen

/// Matrices
template <typename Numeric, int N, int M>
using matrix_t = Eigen::Matrix<Numeric, N, M>;
template <int N, int M>
using rmatrix_t = matrix_t<real_t, N, M>;
typedef rmatrix_t<n_dynamic, n_dynamic> rmat_t;
template <int N, int M>
using cmatrix_t = matrix_t<cplx_t, N, M>;
typedef cmatrix_t<n_dynamic, n_dynamic> cmat_t;

/// Vectors
template <typename Numeric, int N>
using vector_t = Eigen::Matrix<Numeric, N, 1>;
template <typename Numeric, int N>
using vector_view_t = Eigen::Map<const vector_t<Numeric, N>>;
template <int N>
using ivector_t = vector_t<sint_t, N>;
using ivec_t = ivector_t<n_dynamic>;
using ivec3_t = ivector_t<3>;
template <int N>
using rvector_t = vector_t<real_t, N>;
template <int N>
using rvector_view_t = vector_view_t<real_t, N>;
using rvec_t = rvector_t<n_dynamic>;
using rvec3_t = rvector_t<3>;
template <int N>
using cvector_t = vector_t<cplx_t, N>;
using cvec_t = cvector_t<n_dynamic>;
using cvec3_t = cvector_t<3>;

//  hacky?
namespace Eigen {
  template <typename Numeric, int N>
  inline auto begin(vector_t<Numeric, N>& v) noexcept {
    return v.data();
  }
  template <typename Numeric, int N>
  inline auto begin(const vector_t<Numeric, N>& v) noexcept {
    return v.data();
  }
  template <typename Numeric, int N>
  inline auto begin(vector_view_t<Numeric, N>& v) noexcept {
    return v.data();
  }
  template <typename Numeric, int N>
  inline auto begin(const vector_view_t<Numeric, N>& v) noexcept {
    return v.data();
  }
  template <typename Numeric, int N>
  inline auto end(vector_t<Numeric, N>& v) noexcept {
    return v.data() + v.size();
  }
  template <typename Numeric, int N>
  inline auto end(const vector_t<Numeric, N>& v) noexcept {
    return v.data() + v.size();
  }
  template <typename Numeric, int N>
  inline auto end(vector_view_t<Numeric, N>& v) noexcept {
    return v.data() + v.size();
  }
  template <typename Numeric, int N>
  inline auto end(const vector_view_t<Numeric, N>& v) noexcept {
    return v.data() + v.size();
  }
}    // namespace Eigen

namespace stl {
  template <typename T>
  using vector = std::vector<T, alloc::default_init_allocator<T>>;
}

#endif
