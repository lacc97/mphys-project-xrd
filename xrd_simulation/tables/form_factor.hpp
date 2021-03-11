#ifndef XRD_FORM_FACTOR_HPP
#define XRD_FORM_FACTOR_HPP

#include <gslpp/spline.hpp>

#include "types.hpp"

namespace xrd::tables {
  real_t f0(uint_t Z, real_t x) noexcept;
}

#endif    //XRD_FORM_FACTOR_HPP
