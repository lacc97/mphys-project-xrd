#include <numeric>

#include <fmt/format.h>

#include <nlohmann/json.hpp>

#include <math.hpp>
#include <timer.hpp>

#include "basis.hpp"
#include "crystal.hpp"
#include "diffraction.hpp"
#include "lattice.hpp"

using json = nlohmann::json;

/* length:      angstrom
 * temperature: Kelvin */

namespace {
  namespace xray {
    namespace CuKalpha {
      constexpr real_t lambda = 1.5406;
    }
  }    // namespace xray

  namespace crystals {
    const xrd::crystal k_FePt_tetragonal{
      //      xrd::lattice::fcc_tetragonal(2.728, 3.779),   // lit
      //      xrd::lattice::fcc_tetragonal(2.7232835694196504, 3.7065335909902877),   // 0 0 1
      //      xrd::lattice::fcc_tetragonal(2.739400497432593, 3.7106122821630896),    // 0 0 1
      xrd::lattice::fcc_tetragonal(3.2885517977281973, 3.7110427557019965),    // 0 0 1
                                                                               //      xrd::lattice::fcc_tetragonal(3.22229413651982, 3.709930515585622),
      //      xrd::lattice::fcc_tetragonal(2.2899753152966498, 3.710657085181108),    // 0 0 1
      xrd::basis{{26, 55.84, {0, 0, 0}}, {78, 195.08, {0.5, 0.5, 0.5}}},
      230    // https://doi.org/10.1016/0304-8853(94)01586-4
    };

    const xrd::crystal k_KBr{xrd::lattice::fcc(4.740), xrd::basis{{18, 39.098, {0, 0, 0}}, {36, 79.904, {0.5, 0.5, 0.5}}}, 172};

    const xrd::crystal k_KCl{xrd::lattice::fcc(4.514), xrd::basis{{18, 39.098, {0, 0, 0}}, {18, 35.450, {0.5, 0.5, 0.5}}}, 255};
  }    // namespace crystals
}    // namespace

int main(int argc, char** argv) {
  using namespace xray::CuKalpha;

  hr_timer timer{"XRD pattern generation"};

  const xrd::crystal& crystal = crystals::k_FePt_tetragonal;
  const rvec3_t plane = {0, 0, 1};
  const ivector_t<3> size = {40, 40, 40};
  constexpr real_t mspread = math::deg2rad(1);

  constexpr real_t temp = 300;

  xrd::single_plane_diffraction_pattern experiment(crystal, size, mspread, 1000, plane, temp, lambda);
  const rdata_t angles = math::data::linspace(5, 80, 10000);

  timer.start();
  const rdata_t intensities = experiment.generate(angles);
  timer.stop();
  timer.report();

  auto* file = fopen("fept/xrd_pattern.csv", "w");
  for(sint_t ii = 0; ii < intensities.size(); ++ii)
    fmt::print(file, "{} {}\n", angles(ii), intensities(ii));
  fclose(file);
}