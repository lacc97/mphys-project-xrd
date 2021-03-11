#ifndef XRD_SPLINE_HPP
#define XRD_SPLINE_HPP

#include <memory>
#include <span>

#include "interp.hpp"

namespace gsl {
  class spline {
    struct spline_deleter {
      void operator()(void* ptr) const noexcept;
    };
   public:
    constexpr explicit spline(interp_type type) noexcept : m_type{type} {}
    spline(interp_type type, std::span<const double> x, std::span<const double> y) : m_type{type} {
      init(x, y);
    }

    void init(std::span<const double> x, std::span<const double> y);

    double operator()(double x, interp_accel& acc) const noexcept;
    double d(double x, interp_accel& acc) const noexcept;
    double d2(double x, interp_accel& acc) const noexcept;

   private:
    interp_type m_type;
    std::unique_ptr<void, spline_deleter> mp_spline;
  };
}

#endif    //XRD_SPLINE_HPP
