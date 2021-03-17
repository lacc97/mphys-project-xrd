#include "crystal.hpp"

#include <array>

#include <nlohmann/json.hpp>

#include <math.hpp>

#include "tables/form_factor.hpp"

xrd::crystal nlohmann::adl_serializer<xrd::crystal>::from_json(const json& j) {
  if(j.contains("debye_temperature"))
    return xrd::crystal{j.at("lattice").get<xrd::lattice>(), j.at("basis").get<xrd::basis>(), j.at("debye_temperature").get<real_t>()};
  else
    return xrd::crystal{j.at("lattice").get<xrd::lattice>(), j.at("basis").get<xrd::basis>()};
}

void nlohmann::adl_serializer<xrd::crystal>::to_json(nlohmann::json& j, const xrd::crystal& c) {
  j["lattice"] = c.lattice();
  j["basis"] = c.basis();
  if(c.m_DebyeTemperature)
    j["debye_temperature"] = *c.m_DebyeTemperature;
}

real_t xrd::crystal::debye_temperature() const noexcept {
  if(m_DebyeTemperature)
    return *m_DebyeTemperature;
  else
    throw std::runtime_error("unimplemented: xrd::crystal::debye_temperature()");
}

cplx_t xrd::crystal::structure_factor(const rvec3_t& wavevector) const noexcept {
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
    [this, x, &wavevector](const basis::atom& atom) -> structure {
      const rvec3_t r = m_Lattice.r3_vector(atom.r);
      const cplx_t f = tables::f0(atom.f, x);
      return {math::squared_norm(f), f * math::exp(-k_i * wavevector.dot(r))};
    });

  return s / std::sqrt(f);
}