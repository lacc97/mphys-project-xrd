#ifndef XRD_LATTICE_HPP
#define XRD_LATTICE_HPP

#include <utility>

#include <nlohmann/json_fwd.hpp>

#include <constants.hpp>
#include <types.hpp>

namespace xrd {
  class lattice {
   public:
    inline static lattice cubic(const real_t a) noexcept {
      return {a * rvec3_t{1, 0, 0}, a * rvec3_t{0, 1, 0}, a * rvec3_t{0, 0, 1}};
    }

    inline static lattice fcc(const real_t a) noexcept {
      return fcc_tetragonal(a, a);
    }

    inline static lattice fcc_tetragonal(const real_t a, const real_t c) noexcept {
      return {a * rvec3_t{1, 1, 0}.normalized(), a * rvec3_t{1, -1, 0}.normalized(), c * rvec3_t{0, 0, 1}};
    }

    lattice(const rvec3_t& a, const rvec3_t& b, const rvec3_t& c) noexcept : m_Lattice{} {
      m_Lattice.col(0) = a;
      m_Lattice.col(1) = b;
      m_Lattice.col(2) = c;
    }


    [[nodiscard]] real_t cell_volume() const noexcept {
      return a().dot(b().cross(c()));
    }

    [[nodiscard]] inline rvec3_t r3_vector(const rvec3_t& r) const noexcept {
      return m_Lattice * r;
    }

    [[nodiscard]] lattice reciprocal() const noexcept {
      const real_t v = cell_volume();

      return {2 * C_PI * b().cross(c()) / v, 2 * C_PI * c().cross(a()) / v, 2 * C_PI * a().cross(b()) / v};
    }

    [[nodiscard]] inline rvec3_t a() const noexcept {
      return m_Lattice.col(0);
    }
    [[nodiscard]] inline rvec3_t b() const noexcept {
      return m_Lattice.col(1);
    }
    [[nodiscard]] inline rvec3_t c() const noexcept {
      return m_Lattice.col(2);
    }

    [[nodiscard]] inline const rmatrix_t<3, 3>& basis_matrix() const noexcept {
      return m_Lattice;
    }

   private:
    rmatrix_t<3, 3> m_Lattice;
  };
}    // namespace xrd

namespace nlohmann {
  template <>
  struct adl_serializer<xrd::lattice> {
    static xrd::lattice from_json(const json& j);
    static void to_json(json& j, const xrd::lattice& l);
  };
}    // namespace nlohmann

#endif    //XRD_LATTICE_HPP
