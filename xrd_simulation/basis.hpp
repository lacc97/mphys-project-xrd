#ifndef XRD_BASIS_HPP
#define XRD_BASIS_HPP

#include <numeric>
#include <span>

#include <nlohmann/json_fwd.hpp>

#include <types.hpp>

#include "lattice.hpp"

namespace xrd {
  class basis {
   public:
    struct atom {
      uint_t f;
      real_t m;
      rvec3_t r;
    };

    basis(std::initializer_list<atom> atoms) : basis(std::span(atoms)) {}
    basis(std::span<const atom> atoms) : m_Atoms(atoms.begin(), atoms.end()) {}

    [[nodiscard]] real_t total_mass() const noexcept {
      return std::transform_reduce(m_Atoms.begin(), m_Atoms.end(), real_t(0), std::plus<>(), [](const xrd::basis::atom& a) { return a.m; });
    }

    [[nodiscard]] auto begin() const noexcept {
      return m_Atoms.begin();
    }
    [[nodiscard]] auto end() const noexcept {
      return m_Atoms.end();
    }

    [[nodiscard]] auto count() const noexcept {
      return m_Atoms.size();
    }

   private:
    std::vector<atom> m_Atoms;
  };
}    // namespace xrd

namespace nlohmann {
  template <>
  struct adl_serializer<xrd::basis::atom> {
    static void from_json(const json& j, xrd::basis::atom& a);
    static void to_json(json& j, const xrd::basis::atom& a);
  };

  template <>
  struct adl_serializer<xrd::basis> {
    static xrd::basis from_json(const json& j);
    static void to_json(json& j, const xrd::basis& b);
  };
}    // namespace nlohmann

#endif    //XRD_BASIS_HPP
