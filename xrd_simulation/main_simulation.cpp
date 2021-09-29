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

namespace {
  real_t background_profile(real_t theta_deg) noexcept {
    constexpr real_t A = 307.369;
    constexpr real_t Lambda = 0.0400355;
    constexpr real_t Sigma = 0.0281618;
    constexpr real_t C = 16.4488;

    constexpr real_t a_1 = 282640;
    constexpr real_t mu_1 = 42.9251;
    constexpr real_t sigma_1 = 0.0119425;
    constexpr real_t gamma_1 = 0.00534363;

    constexpr real_t a_2 = 83695.3;
    constexpr real_t mu_2 = 94.0501;
    constexpr real_t sigma_2 = 0.0104494;
    constexpr real_t gamma_2 = 0.0208413;

    const real_t x = math::deg2rad(theta_deg);
    const real_t bg_lin = A*std::exp(-Lambda*x - (Sigma*x)*(Sigma*x)) + C;
    const real_t bg_1 = a_1*math::stats::voigt(x - mu_1, sigma_1, gamma_1);
    const real_t bg_2 = a_2*math::stats::voigt(x - mu_2, sigma_2, gamma_2);

    const real_t bg = bg_lin + bg_1 + bg_2;
    return bg;
  }
}    // namespace

int main(int argc, char** argv) {
#if 1
  if(argc != 2)
    throw std::runtime_error("need to provide .json file as first argument");

  json config = io::load_json(argv[1]);

  real_t global_factor;
  uint_t mosaic_samples;
  rdata_t angles;
  bool with_bg;
  {
    const auto& c_env = config.at("computational_environment");

    std::array<real_t, 2> angle_interval;
    c_env.at("angle_interval").get_to(angle_interval);
    angles = math::data::linspace(angle_interval[0], angle_interval[1], c_env.at("angle_samples").get<uint_t>());

    c_env.at("mosaic_samples").get_to(mosaic_samples);

    global_factor = c_env.contains("global_factor") ? c_env.at("global_factor").get<real_t>() : 1;

    with_bg = c_env.contains("with_bg") ? c_env.at("with_bg").get<bool>() : false;
  }

  real_t wavelength, temperature, slit_angle;
  {
    const auto& p_env = config.at("physical_environment");

    p_env.at("wavelength").get_to(wavelength);
    p_env.at("temperature").get_to(temperature);
    slit_angle = math::deg2rad(p_env.at("receiving_slit_angle").get<real_t>());
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

        xrd::single_plane_diffraction_pattern experiment(crystal, crystallite_size, mosaic_spread, mosaic_samples, plane, temperature, wavelength, slit_angle);

        rdata_t e_pat = experiment.generate(angles);
        {
          ds::dataset_2d_view dset = ds::dataset_2d_view(angles, e_pat, ds::no_validation);
          fmt::print("  {0}: {{{1}}}\n", plane, fmt::join(dset.find_peaks(0.1), ", "));
        }
        xrd_pattern += multiplicity * e_pat;
      }
    }
  }

  std::string output_path = config.at("output_path").get<std::string>();

  if(with_bg)
    xrd_pattern = (global_factor*xrd_pattern) + angles.unaryExpr(&background_profile);
  else
    xrd_pattern *= global_factor;
  io::write_csv(fmt::format("{0}.csv", output_path), std::tie(angles, xrd_pattern));
#else
  ds::dataset_2d real_pattern(io::load_csv("fept/FePt_XRD.csv"));
  rdata_t bg_removed = real_pattern.x().unaryExpr(&background_profile);
  io::write_csv("fept/new_stuff/minus_bg.csv", real_pattern.x(), bg_removed);
#endif
}