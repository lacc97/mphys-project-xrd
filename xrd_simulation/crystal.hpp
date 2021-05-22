#ifndef XRD_CRYSTAL_HPP
#define XRD_CRYSTAL_HPP

#include <optional>

#include <nlohmann/json_fwd.hpp>

#include <math.hpp>

#include "basis.hpp"
#include "lattice.hpp"

#include "tables/form_factor.hpp"

namespace xrd {
  class crystal {
    friend struct ::nlohmann::adl_serializer<crystal>;

   public:
    crystal(xrd::lattice l, xrd::basis b) : m_Lattice{std::move(l)}, m_Basis{std::move(b)} {}
    crystal(xrd::lattice l, xrd::basis b, real_t debye) : m_Lattice{std::move(l)}, m_Basis{std::move(b)}, m_DebyeTemperature{debye} {}

    [[nodiscard]] real_t debye_temperature() const noexcept;
    [[nodiscard]] real_t mass_density() const noexcept {
      return m_Basis.total_mass() / m_Lattice.cell_volume();
    }
    [[nodiscard]] real_t number_density() const noexcept {
      return m_Basis.count() / m_Lattice.cell_volume();
    }
    [[nodiscard]] inline cplx_t structure_factor(const rvec3_t& wavevector) const noexcept {
      return structure_factor(wavevector, [](const basis::atom& atom) { return 1; });
    }

    template <typename F>
    [[nodiscard]] cplx_t structure_factor(const rvec3_t& wavevector, F&& ff_mod) const noexcept {
      using namespace std::complex_literals;

      struct structure {
        real_t f;
        cplx_t s;
      };

      auto x = wavevector.norm()/(4*C_PI);
      auto [f, s] = std::transform_reduce(
        m_Basis.begin(), m_Basis.end(), structure{},
        [](const structure& p1, const structure& p2) -> structure {
          return {p1.f + p2.f, p1.s + p2.s};
        },
        [this, x, &wavevector, &ff_mod](const basis::atom& atom) -> structure {
          const rvec3_t r = m_Lattice.r3_vector(atom.r);
          const cplx_t f = tables::f0(atom.f, x) * ff_mod(atom);
          return {math::squared_norm(f), f * math::exp(-k_i * wavevector.dot(r))};
        });

      return s / std::sqrt(f);
    }

    [[nodiscard]] inline const xrd::basis& basis() const noexcept {
      return m_Basis;
    }
    [[nodiscard]] inline const xrd::lattice& lattice() const noexcept {
      return m_Lattice;
    }

   private:
    xrd::lattice m_Lattice;
    xrd::basis m_Basis;
    std::optional<real_t> m_DebyeTemperature;
  };
}    // namespace xrd

namespace nlohmann {
  template <>
  struct adl_serializer<xrd::crystal> {
    static xrd::crystal from_json(const json& j);
    static void to_json(json& j, const xrd::crystal& c);
  };
}    // namespace nlohmann

#endif    //XRD_CRYSTAL_HPP
