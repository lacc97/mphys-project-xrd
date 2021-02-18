#include "lattice.hpp"

#include <fmt/format.h>

#include <nlohmann/json.hpp>

#include <types_json.hpp>

xrd::lattice nlohmann::adl_serializer<xrd::lattice>::from_json(const nlohmann::json& j) {
  if(j.contains("type")) {
    auto type = j.at("type").get<std::string>();
    if(type == "cubic")
      return xrd::lattice::cubic(j.at("a").get<real_t>());
    else if(type == "fcc")
      return xrd::lattice::fcc(j.at("a").get<real_t>());
    else if(type == "fcc_tetragonal")
      return xrd::lattice::fcc_tetragonal(j.at("a").get<real_t>(), j.at("c").get<real_t>());
    else
      throw std::runtime_error(fmt::format("unrecognized lattice type: {}", type));
  } else {
    return xrd::lattice(j.at("a").get<rvec3_t>(), j.at("b").get<rvec3_t>(), j.at("c").get<rvec3_t>());
  }
}

void nlohmann::adl_serializer<xrd::lattice>::to_json(nlohmann::json& j, const xrd::lattice& l) {
  j["a"] = l.a();
  j["b"] = l.b();
  j["c"] = l.c();
}
