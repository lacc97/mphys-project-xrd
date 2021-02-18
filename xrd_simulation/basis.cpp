#include "basis.hpp"

#include <nlohmann/json.hpp>

#include <types_json.hpp>

void nlohmann::adl_serializer<xrd::basis::atom>::from_json(const nlohmann::json& j, xrd::basis::atom& a) {
  j.at("form").get_to(a.f);
  j.at("mass").get_to(a.m);
  j.at("position").get_to(a.r);
}

void nlohmann::adl_serializer<xrd::basis::atom>::to_json(nlohmann::json& j, const xrd::basis::atom& a) {
  j["form"] = a.f;
  j["mass"] = a.m;
  j["position"] = a.r;
}

xrd::basis nlohmann::adl_serializer<xrd::basis>::from_json(const json& j) {
  stl::vector<xrd::basis::atom> atoms(j.size());

  for(uint_t ii = 0; ii < atoms.size(); ++ii)
    j[ii].get_to(atoms[ii]);

  return xrd::basis(atoms);
}

void nlohmann::adl_serializer<xrd::basis>::to_json(nlohmann::json& j, const xrd::basis& b) {
  for(const auto& a : b)
    j.push_back(a);
}
