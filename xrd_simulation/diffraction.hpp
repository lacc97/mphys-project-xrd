#ifndef XRD_DIFFRACTION_HPP
#define XRD_DIFFRACTION_HPP

#include <types.hpp>

#include "basis.hpp"
#include "crystal.hpp"
#include "lattice.hpp"

namespace xrd {
  class single_plane_diffraction_pattern {
   public:
    single_plane_diffraction_pattern(xrd::crystal c, ivector_t<3> c_size, real_t m_spread, uint_t m_samples, rvec3_t plane, real_t temp, real_t wavelength, real_t rec_slit)
        : m_Crystal{std::move(c)}, m_ReciprocalLattice{m_Crystal.lattice().reciprocal()}, m_CrystalliteSize{std::move(c_size)}, m_MosaicSpread{m_spread},
          m_MosaicSamples{m_samples}, m_Plane{std::move(plane)}, m_Temperature{temp}, m_XrayWavelength{wavelength}, m_ReceivingSollerSlitAngle{rec_slit} {}

    [[nodiscard]] inline rdata_t generate(const rdata_t& angles) const {
      rdata_t intensities(angles.size());
      generate(angles, intensities);
      return intensities;
    }
    void generate(const rdata_t& angles, rdata_t& intensities) const;

   private:
    [[nodiscard]] rmatrix_t<3, n_dynamic> generate_random_plane_vectors() const;


    xrd::crystal m_Crystal;
    xrd::lattice m_ReciprocalLattice;
    ivector_t<3> m_CrystalliteSize;
    real_t m_MosaicSpread;
    uint_t m_MosaicSamples;

    rvec3_t m_Plane;

    real_t m_Temperature;
    real_t m_XrayWavelength;
    real_t m_ReceivingSollerSlitAngle;
  };

  inline real_t scherrer_factor(const lattice& latt, const ivector_t<3>& sizes, const rvec3_t& wavevector) {
    auto fn_xi = [](sint_t N, real_t x) noexcept -> real_t {
      const real_t sin_x = std::sin(x);
      const real_t f = (sin_x == 0) ? (N) : (std::sin(N * x) / sin_x);
      return f * f;
    };

    auto fn_fac = [&latt, &sizes, &wavevector, &fn_xi](size_t ii) noexcept -> real_t {
      const sint_t N = sizes(ii);
      const rvec3_t a = latt.basis_matrix().col(ii);

      return fn_xi(N, wavevector.dot(a) / 2) / (real_t(N) * N);
    };

    return fn_fac(0) * fn_fac(1) * fn_fac(2);
  }

  inline real_t lorentz_factor(real_t angle) noexcept {
    real_t sin = std::sin(angle);
    return 1 / (4 * sin * sin * std::cos(angle));
  }

  inline real_t polarization_factor(real_t angle) noexcept {
    real_t cos = std::cos(2 * angle);
    return (1 + cos * cos) / 2;
  }

  real_t powder_grain_distribution_factor(real_t a, real_t sigma) noexcept;
}    // namespace xrd

#endif    //XRD_DIFFRACTION_HPP
