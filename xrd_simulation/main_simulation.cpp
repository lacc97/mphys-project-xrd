#include <numeric>

#include <fmt/format.h>

#include <nlohmann/json.hpp>

#include <data/dataset_2d.hpp>
#include <io.hpp>
#include <math.hpp>
#include <timer.hpp>
#include <types_format.hpp>
#include <types_json.hpp>

#include "basis.hpp"
#include "crystal.hpp"
#include "diffraction.hpp"
#include "lattice.hpp"

using json = nlohmann::json;

template <>
struct fmt::formatter<ds::dataset_2d_view::point> {
  auto parse(format_parse_context& ctx) {
    return ctx.begin();
  }

  auto format(const ds::dataset_2d_view::point& p, format_context& ctx) {
    return format_to(ctx.out(), "({0}, {1})", p.x, p.y);
  }
};

/* length:      angstrom
 * temperature: Kelvin */

namespace {
  namespace xray {
    namespace CuKalpha {
      constexpr real_t lambda = 1.5406;
    }
  }    // namespace xray

  namespace crystals {
    const xrd::crystal k_KBr{xrd::lattice::fcc(4.740), xrd::basis{{18, 39.098, {0, 0, 0}}, {36, 79.904, {0.5, 0.5, 0.5}}}, 172};

    const xrd::crystal k_KCl{xrd::lattice::fcc(4.514), xrd::basis{{18, 39.098, {0, 0, 0}}, {18, 35.450, {0.5, 0.5, 0.5}}}, 255};
  }    // namespace crystals
}    // namespace

int main(int argc, char** argv) {
  if(argc != 2)
    throw std::runtime_error("need to provide .json file as first argument");

  json config = io::load_json(argv[1]);

  uint_t mosaic_samples;
  rdata_t angles;
  {
    const auto& c_env = config.at("computational_environment");

    std::array<real_t, 2> angle_interval;
    c_env.at("angle_interval").get_to(angle_interval);
    angles = math::data::linspace(angle_interval[0], angle_interval[1], c_env.at("angle_samples").get<uint_t>());

    c_env.at("mosaic_samples").get_to(mosaic_samples);
  }

  real_t wavelength, temperature;
  {
    const auto& p_env = config.at("physical_environment");

    p_env.at("wavelength").get_to(wavelength);
    p_env.at("temperature").get_to(temperature);
  }

  rdata_t xrd_pattern = rdata_t::Zero(angles.size());
  for(const auto& c : config.at("crystals")) {
    xrd::crystal crystal = c;

    std::string name = c.contains("name") ? c.at("name") : "";

    ivec3_t crystallite_size = c.at("crystallite_size").get<ivec3_t>();
    real_t mosaic_spread = math::deg2rad(c.at("mosaic_spread").get<real_t>());

    const auto& p = c.at("patterns");
    fmt::print("Crystal: {0}\n", name);
    for(const auto& p_config : p) {
      real_t multiplicity = p_config.contains("multiplicity") ? p_config.at("multiplicity").get<real_t>() : 1;
      if(multiplicity != 0) {
        rvec3_t plane = p_config.at("plane").get<rvec3_t>();

        xrd::single_plane_diffraction_pattern experiment(crystal, crystallite_size, mosaic_spread, mosaic_samples, plane, temperature, wavelength);

        rdata_t e_pat = experiment.generate(angles);
        e_pat *= multiplicity;
        xrd_pattern += e_pat;
        {
          ds::dataset_2d_view dset = ds::dataset_2d_view(angles, e_pat, ds::no_validation);
          fmt::print("  {0}: {{{1}}}\n", plane, fmt::join(dset.find_peaks(), ", "));
        }
      }
    }
  }

  std::string output_path = config.at("output_path").get<std::string>();

  io::write_csv(fmt::format("{0}.csv", output_path), std::tie(angles, xrd_pattern));
}