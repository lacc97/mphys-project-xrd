#include "diffraction.hpp"

#include <gslpp/spline.hpp>

#include <constants.hpp>
#include <math.hpp>

void xrd::single_plane_diffraction_pattern::generate(const rdata_t& angles, rdata_t& intensities) const {
  rmatrix_t<3, Eigen::Dynamic> delta_k_vectors;
  if(m_MosaicSamples == 0 || m_MosaicSpread == 0)
    delta_k_vectors = m_ReciprocalLattice.r3_vector(m_Plane).normalized();
  else
    delta_k_vectors = generate_random_plane_vectors();
  delta_k_vectors *= 2 * (2 * C_PI / m_XrayWavelength);

  real_t debye = m_Crystal.debye_temperature();
  auto fn_diff_intensity = [&, x = debye / m_Temperature](real_t theta) noexcept -> real_t {
    const real_t sin_theta = std::sin(theta);
    const real_t cos_theta = std::cos(theta);
    const real_t cos_2theta = std::cos(2 * theta);

    auto fn_debye_waller = [this, debye, x, x_2 = x*x, sin_theta](const xrd::basis::atom& atom) noexcept -> real_t {
      /* This is only correct for cubic crystals */
      constexpr real_t C1 = 1.915045e3;    // h^2 / k_b considering units of length are angstroms, units of mass are daltons

      const real_t C2 = 1 + (x_2 / 36) - (x_2 * x_2 / 3600);    // Taylor expansion around 0 of [phi(x) + x/4]
      return std::exp(-2 * ((3 * C1) / (2 * C_PI * C_PI * atom.m * x * debye)) * C2 * 2 * C_PI * sin_theta * sin_theta / (m_XrayWavelength*m_XrayWavelength));
    };

    const real_t f_lorentz = 1 / (4 * sin_theta * sin_theta * cos_theta);
    const real_t f_polarization = (1 + cos_2theta * cos_2theta) / 2;

    const real_t f_temperature = std::transform_reduce(m_Crystal.basis().begin(), m_Crystal.basis().end(), (real_t) 1, std::multiplies<>(), fn_debye_waller);

    real_t intensity = 0;
    for(sint_t ii = 0; ii < delta_k_vectors.cols(); ++ii) {
      const rvec3_t delta_k = delta_k_vectors.col(ii) * sin_theta;

      const real_t f_geometry = xrd::scherrer_factor(m_Crystal.lattice(), m_CrystalliteSize, delta_k);

      const real_t factors = f_lorentz * f_polarization * f_geometry * f_temperature;

      intensity += math::squared_norm(m_Crystal.structure_factor(delta_k)) * factors;
    }
    return intensity / delta_k_vectors.cols();
  };

  intensities.resize(angles.size());

#pragma omp parallel for default(none) shared(fn_diff_intensity, intensities, angles)
  for(sint_t ii = 0; ii < intensities.size(); ++ii)
    intensities(ii) = fn_diff_intensity(math::deg2rad(angles(ii)));
}

rmatrix_t<3, n_dynamic> xrd::single_plane_diffraction_pattern::generate_random_plane_vectors() const {
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

  return vectors;
}
