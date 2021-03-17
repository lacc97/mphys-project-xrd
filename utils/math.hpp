#ifndef XRD_MATH_HPP
#define XRD_MATH_HPP

#include <random>

#include <fmt/format.h>

#include <gsl/gsl_math.h>

#include "constants.hpp"
#include "data/dataset_2d.hpp"
#include "types.hpp"

namespace math {
  template <typename T>
  concept number_type = requires (T a, T b) {
    {a + b} -> std::same_as<T>;
    {a - b} -> std::same_as<T>;
    {a * b} -> std::same_as<T>;
    {a / b} -> std::same_as<T>;
  };

  namespace data {
    namespace details {
      void linspace(real_t start, real_t stop, std::span<real_t> out);

      template <typename F>
      void table(real_t start, real_t stop, std::span<real_t> x, std::span<real_t> y, F&& f) {
        if(x.size() != y.size())
          throw std::invalid_argument(fmt::format("shape mismatch: x size ({}) != y size ({})", x.size(), y.size()));
        if(x.size() <= 1)
          throw std::invalid_argument(fmt::format("invalid count ({}): must be at least 2", x.size()));
        if(stop <= start)
          throw std::invalid_argument(fmt::format("invalid endpoints: start ({}) must be smaller than stop ({})", start, stop));

        const real_t h = (stop - start) / static_cast<real_t>(x.size() - 1);

        for(size_t ii = 0; ii < x.size() - 1; ++ii) {
          x[ii] = start + ii * h;
          y[ii] = f(x[ii]);
        }
        x.back() = stop;
        y.back() = f(x.back());
      }
    }    // namespace details

    inline rdata_t linspace(real_t start, real_t stop, uint_t count) {
      rdata_t space(count);
      details::linspace(start, stop, space);
      return space;
    }

    template <typename F>
    ds::dataset_2d table(real_t start, real_t stop, uint_t count, F&& f) {
      stl::vector<real_t> x(count), y(count);
      details::table(start, stop, x, y, std::forward<F>(f));
      return ds::dataset_2d(std::move(x), std::move(y));
    }
  }    // namespace data

  namespace linalg {
    rmatrix_t<3, 3> generate_orthonormal_basis_from_vector(const rvec3_t& z_axis) noexcept;
  }

  namespace rand {
    extern thread_local std::mt19937_64 tl_Generator;

    real_t unit();
  }    // namespace rand

  template <typename T>
  inline constexpr std::enable_if_t<std::is_arithmetic_v<T>, int> signum(T x) noexcept {
    if constexpr(std::is_signed<T>())
      return (T(0) < x) - (x < T(0));
    else
      return (T(0) < x);
  }

  inline constexpr real_t deg2rad(real_t deg) noexcept {
    return deg * C_PI / real_t(180.0);
  }
  inline constexpr real_t rad2deg(real_t rad) noexcept {
    return rad * real_t(180.0) / C_PI;
  }

  inline constexpr cplx_t conj(cplx_t z) noexcept {
    return z.conj();
  }
  inline constexpr real_t conj(real_t r) noexcept {
    return r;
  }

  inline cplx_t exp(cplx_t z) noexcept {
    return gsl_complex_exp(static_cast<gsl_complex>(z));
  }
  inline real_t exp(real_t r) noexcept {
    return std::exp(r);
  }

  inline constexpr real_t norm(cplx_t z) noexcept {
    return z.norm();
  }
  inline constexpr real_t norm(real_t r) noexcept {
    return r;
  }

  inline constexpr real_t squared_norm(cplx_t z) noexcept {
    return z.squared_norm();
  }
  inline constexpr real_t squared_norm(real_t r) noexcept {
    return r*r;
  }

  template <number_type T>
  inline constexpr T sqr(T r) noexcept {
    return r * r;
  }
}    // namespace math

#include "math/peak_finder.hpp"

#endif    //XRD_MATH_HPP