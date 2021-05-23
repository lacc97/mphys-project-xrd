#ifndef XRD_INTEGRATION_HPP
#define XRD_INTEGRATION_HPP

#include <cstddef>

#include <memory>
#include <type_traits>

namespace gsl::integration {
  template <typename F>
  concept function = std::is_invocable_r_v<double, F, double>;

  enum class rule { e_Gauss_15, e_Gauss_21, e_Gauss_31, e_Gauss_41, e_Gauss_51, e_Gauss_61 };

  struct result {
    operator double() const noexcept {
      return res;
    }

    double res, abserr;
  };

  namespace detail {
    struct function {
      double (*function)(double x, void* params);
      const void* params;
    };

    template <integration::function F>
    inline function to_function(F&& f) noexcept {
      return {+[](double x, void* params) -> double { return (*static_cast<const F*>(params))(x); }, static_cast<const void*>(std::addressof(f))};
    }

    result qag_adaptive(function f, double a, double b, rule r, double epsabs, double epsrel, size_t limit);
  }    // namespace detail


  template <typename F> requires function<F>
  inline result qag_adaptive(F&& f, double a, double b, rule r, double epsabs, double epsrel, size_t limit = 0) {
    return detail::qag_adaptive(detail::to_function(std::forward<F>(f)), a, b, r, epsabs, epsrel, limit);
  }
}    // namespace gsl::integration

#endif    //XRD_INTEGRATION_HPP
