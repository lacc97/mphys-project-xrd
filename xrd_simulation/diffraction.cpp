#include "diffraction.hpp"

#include <gslpp/integration.hpp>
#include <gslpp/spline.hpp>

#include <constants.hpp>
#include <math.hpp>

namespace {
  real_t temp_dimensionless_phi(real_t x) {
    real_t d = gsl::integration::qag_adaptive([x](double e) -> double { return e / (std::exp(e) - 1); }, 0, x, gsl::integration::rule::e_Gauss_41, 0, 1e-5);
    return d / x;
  }

  real_t temp_v2(real_t debye, real_t T) {
    constexpr real_t C1 = 145.526; /*3*hb/k_b in K.Da.A^2 */

    const real_t x = debye / T;
    return C1 * ((temp_dimensionless_phi(x) / x) + 0.25) / debye;
  }
}    // namespace

struct xrd::single_plane_diffraction_pattern::workspace {
  explicit workspace(rmatrix_t<3, Eigen::Dynamic> mosaics, real_t s2, real_t debye, real_t T)
      : mosaic_planes{std::move(mosaics)}, tan_s2{std::tan(s2)}, x{debye / T}, x_2{x}, phi_x{temp_dimensionless_phi(x)}, c2{phi_x + x / 4}, v2{temp_v2(debye,
                                                                                                                                                       T)} {}

  const rmatrix_t<3, Eigen::Dynamic> mosaic_planes;

  const real_t tan_s2;

  const real_t x;
  const real_t x_2;
  const real_t phi_x;
  const real_t c2;

  const real_t v2;
};

real_t xrd::single_plane_diffraction_pattern::calculate_intensity_internal(const xrd::single_plane_diffraction_pattern::workspace& w, real_t theta) const {
  const real_t sin_theta = std::sin(theta);
  const real_t csc_theta = 1 / sin_theta;
  const real_t cos_theta = std::cos(theta);
  const real_t sin_2theta = 2 * sin_theta * cos_theta;
  const real_t cos_2theta = std::cos(2 * theta);

  auto fn_f = [this, &w, sin_theta](const xrd::basis::atom& a) -> real_t {
    return math::exp(-8 * C_PI * C_PI * (w.v2 / a.m) * (sin_theta / m_XrayWavelength) * (sin_theta / m_XrayWavelength));
  };

  const real_t f_abs = (1 - std::exp(-2 * m_AbsorptionUT / sin_theta));

  const real_t f_lorentz = 1 / (2 * sin_theta * sin_2theta);
  const real_t f_polarization = (1 + cos_2theta * cos_2theta) / 2;

  real_t intensity = 0;
  for(sint_t ii = 0; ii < w.mosaic_planes.cols(); ++ii) {
    const rvec3_t delta_k = w.mosaic_planes.col(ii) * sin_theta;

    const real_t f_geometry = xrd::scherrer_factor(m_Crystal.lattice(), m_CrystalliteSize, delta_k);

    const real_t factors = f_lorentz * f_polarization * f_geometry * f_abs;

    intensity += math::squared_norm(m_Crystal.structure_factor(delta_k, fn_f)) * factors;
  }
  return intensity / w.mosaic_planes.cols();
}

real_t xrd::single_plane_diffraction_pattern::calculate_intensity_with_mosaic(rmatrix_t<3, n_dynamic> mosaic_planes, real_t theta) const {
  workspace w{std::move(mosaic_planes), m_ReceivingSollerSlitAngle, m_Crystal.debye_temperature(), m_Temperature};
  return calculate_intensity_internal(w, theta);
}

void xrd::single_plane_diffraction_pattern::generate(const rdata_t& angles, rdata_t& intensities) const {
  intensities.resize(angles.size());

  workspace w{generate_random_scattering_vectors(), m_ReceivingSollerSlitAngle, m_Crystal.debye_temperature(), m_Temperature};
#pragma omp parallel for default(none) shared(w, angles, intensities)
  for(sint_t ii = 0; ii < intensities.size(); ++ii)
    intensities(ii) = calculate_intensity_internal(w, math::deg2rad(angles(ii)));
}

rmatrix_t<3, n_dynamic> xrd::single_plane_diffraction_pattern::generate_random_scattering_vectors() const {
  real_t magnitude = (2 * (2 * C_PI / m_XrayWavelength));
  if(m_MosaicSamples == 0 || m_MosaicSpread == 0) {
    return magnitude * m_ReciprocalLattice.r3_vector(m_Plane).normalized();
  } else {
    rmatrix_t<3, n_dynamic> vectors(3, m_MosaicSamples);

    const rvec3_t g = m_ReciprocalLattice.r3_vector(m_Plane);
    rmatrix_t<3, 3> basis = math::linalg::generate_orthonormal_basis_from_vector(g.normalized());

    auto fn_rotate = [&basis](real_t phi, real_t theta) noexcept -> rvec3_t {
      const real_t sin_theta = std::sin(theta);

      const rvec3_t v{std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, std::cos(theta)};

      return basis * v;
    };

    std::uniform_real_distribution<real_t> phi_dist{0, 2 * C_PI};
    std::normal_distribution<real_t> theta_dist{0, m_MosaicSpread};
    for(sint_t ii = 0; ii < vectors.cols(); ++ii) {
      real_t phi = phi_dist(math::rand::tl_Generator), theta = theta_dist(math::rand::tl_Generator);
      vectors.col(ii) = fn_rotate(phi, std::abs(theta));
    }

    return magnitude * vectors;
  }
}

real_t xrd::powder_grain_distribution_factor(real_t a, real_t sigma) noexcept {
  //  return 0.5;
  return a;
  //  return math::erf(a/(C_SQRT2*sigma));
}
