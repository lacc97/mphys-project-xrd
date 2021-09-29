#ifndef XRD_MATH_HPP
#define XRD_MATH_HPP

#include <random>

#include <fmt/format.h>

#include <gslpp/integration.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

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
      void logspace(real_t start, real_t stop, std::span<real_t> out);

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

    inline rdata_t logspace(real_t start, real_t stop, uint_t count) {
      rdata_t space(count);
      details::logspace(start, stop, space);
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

  namespace stats {
    inline real_t voigt(real_t x, real_t sigma, real_t gamma) noexcept {
      return gsl::integration::qagi_double([x, sigma, gamma](double xp) -> double {
        return gsl_ran_gaussian_pdf(xp, sigma)*gsl_ran_cauchy_pdf(x - xp, gamma);
      }, 0, 1e-5);
    }
  }

  template <typename T> requires std::is_arithmetic_v<T>
  inline constexpr int signum(T x) noexcept {
    if constexpr(std::is_signed_v<T>)
      return (T(0) < x ? 1 : 0) - (x < T(0) ? 1 : 0);
    else
      return (T(0) < x ? 1 : 0);
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

  inline real_t erf(real_t r) noexcept {
    return gsl_sf_erf(r);
  }

  inline cplx_t exp(cplx_t z) noexcept {
    return gsl_complex_exp(static_cast<gsl_complex>(z));
  }
  inline real_t exp(real_t r) noexcept {
    return std::exp(r);
  }
  inline real_t expm1(real_t r) noexcept {
    return std::expm1(r);
  }
  inline cplx_t expm1(cplx_t z) noexcept {
    cplx_t t = gsl_complex_tanh(static_cast<gsl_complex>(z/2));
    return 2*t/(1-t);
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
