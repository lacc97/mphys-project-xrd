#ifndef XRD_INTERP_P_HPP
#define XRD_INTERP_P_HPP

#include <gsl/gsl_interp.h>

#include "interp.hpp"

namespace gsl::details {
  inline const gsl_interp_type* get_type(interp_type type) noexcept {
    switch(type) {
      case interp_type::linear:
        return gsl_interp_linear;
      case interp_type::polynomial:
        return gsl_interp_polynomial;
      case interp_type::cspline:
        return gsl_interp_cspline;
      case interp_type::cspline_periodic:
        return gsl_interp_cspline_periodic;
      case interp_type::akima:
        return gsl_interp_akima;
      case interp_type::akima_periodic:
        return gsl_interp_akima_periodic;
      case interp_type::steffen:
        return gsl_interp_steffen;
      default:
        return nullptr;
    }
  }
}

#endif    //XRD_INTERP_P_HPP
