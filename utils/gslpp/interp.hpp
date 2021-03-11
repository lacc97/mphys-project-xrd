#ifndef XRD_INTERP_HPP
#define XRD_INTERP_HPP

#include <memory>
#include <span>

namespace gsl {
  class spline;

  enum class interp_type {
    linear,
    polynomial,
    cspline,
    cspline_periodic,
    akima,
    akima_periodic,
    steffen
  };

  class interp_accel {
    struct interp_accel_deleter {
      void operator()(void* ptr) const noexcept;
    };

    friend class spline;

   public:
    interp_accel();

    size_t find(std::span<double> array, double x) noexcept;
    void reset();

   private:
    std::unique_ptr<void, interp_accel_deleter> mp_accel;
  };
}

#endif    //XRD_INTERP_HPP
