#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cmath>

#include "types.hpp"

#ifdef M_PIl
constexpr real_t C_PI = M_PIl;
#else
constexpr real_t C_PI = 3.141592653589793238462643383279502884L;
#endif

#if defined(M_SQRT2l)
constexpr real_t C_SQRT2 = M_SQRT2l;
#else
constexpr real_t C_SQRT2 = 1.414213562373095048801688724209698078L;
#endif

namespace SI {
  constexpr real_t C_EV = 1.60217662e-19L;

  constexpr real_t C_C = 299792458;

  constexpr real_t C_PLANCK = 6.62607015e-34L;
  constexpr real_t C_HBAR = C_PLANCK / (2 * C_PI);

  constexpr real_t C_BIG_G = 6.674e-11L;

  constexpr real_t C_ELECTRON_MASS = 5.485799090e-4L * 1.6605390e-27L;
  constexpr real_t C_PROTON_MASS = 1.007276466L * 1.6605390e-27L;
  constexpr real_t C_NEUTRON_MASS = 1.008664915L * 1.6605390e-27L;
  constexpr real_t C_NUCLEON_MASS = (C_PROTON_MASS + C_NEUTRON_MASS) / 2;

  constexpr real_t C_SOLAR_MASS = 1.9884e30L;
  constexpr real_t C_SOLAR_RADIUS = 6.95700e8L;

  //     namespace astro {
  //         constexpr real_t C_SOLAR_MASS =
  //     };
}    // namespace SI

#endif
