#include "crystal.hpp"

#include <array>

#include <nlohmann/json.hpp>

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
