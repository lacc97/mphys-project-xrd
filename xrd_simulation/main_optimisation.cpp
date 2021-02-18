#include <numeric>

#include <fmt/format.h>

#include <data/dataset_2d.hpp>
#include <io.hpp>
#include <math.hpp>
#include <optimisation/simulated_annealing.hpp>
#include <timer.hpp>

#include "basis.hpp"
#include "crystal.hpp"
#include "diffraction.hpp"
#include "lattice.hpp"

/* length:      angstrom
 * temperature: Kelvin */

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
  namespace xray {
    namespace CuKalpha {
      constexpr real_t lambda = 1.5406;
    }
  }    // namespace xray

  class xrd_annealing_simulation {
   public:
    using solution_type = std::array<real_t, 2>;


    xrd_annealing_simulation(ivector_t<3> size, real_t mspread, real_t temp, real_t lambda, rdata_t angles)
        : m_Size{std::move(size)}, m_MosaicSpread{mspread}, m_Temperature{temp}, m_Wavelength{lambda}, m_Angles{std::move(angles)} {
      ds::dataset_2d pattern(io::load_csv("fept/AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta_signal.txt"));

      m_Peak_001 = pattern.get(22.5, 27.5).find_peaks()[0].x;
      m_Peak_110 = pattern.get(30, 35).find_peaks()[0].x;
      m_Peak_111 = pattern.get(40, 42).find_peaks()[0].x;
      m_Peak_200 = pattern.get(46, 48).find_peaks()[0].x;
      m_Peak_002 = pattern.get(47.5, 50).find_peaks()[0].x;

      fmt::print("FePt(001): {{{0}, {5}}}\nFePt(110): {{{1}}}\nFePt(111): {{{2}}}\nMgO(001): {{{3}}}\nFePt(200): {{{4}}}\n", pattern.get(22.5, 27.5).find_peaks()[0], pattern.get(30, 35).find_peaks()[0], pattern.get(40, 42).find_peaks()[0], pattern.get(42, 45).find_peaks()[0], pattern.get(46, 48).find_peaks()[0], pattern.get(47.5, 50).find_peaks()[0]);
    }

    [[nodiscard]] inline solution_type initial_solution() const noexcept {
      //            real_t r0 = (2 * math::rand::unit() - 1) * 0.25, r1 = (2 * math::rand::unit() - 1) * 0.25;
      real_t r0 = (2 * math::rand::unit() - 1) * 0.01, r1 = (2 * math::rand::unit() - 1) * 0.01;
      return {{(1 + r0) * 3.85, (1 + r1) * 3.71}};
      //      return {{(1 + r0) * 2.728, (1 + r1) * 3.779}};
    }

    [[nodiscard]] real_t energy(const solution_type& s) const noexcept {
      const xrd::crystal strained_crystal{xrd::lattice::fcc_tetragonal(s[0], s[1]), xrd::basis{{26, 55.84, {0, 0, 0}}, {78, 195.08, {0.5, 0.5, 0.5}}}, 230};

      auto peaks_001 = find_peak_positions_for_plane(strained_crystal, {0, 0, 1});
      if(peaks_001.size() < 2)
        return std::numeric_limits<real_t>::max();
      real_t e_001 = math::sqr(m_Peak_001 - 2 * peaks_001[0]) + math::sqr(m_Peak_002 - 2 * peaks_001[1]);

      return e_001 + secondary_peaks_energy(strained_crystal);
    }

    [[nodiscard]] inline solution_type random_neighbour(const solution_type& s) const noexcept {
      //            real_t r0 = (2 * math::rand::unit() - 1) * 0.05, r1 = (2 * math::rand::unit() - 1) * 0.05;
      real_t r0 = (2 * math::rand::unit() - 1) * 0.0005, r1 = (2 * math::rand::unit() - 1) * 0.0005;
      return {{(1 + r0) * s[0], (1 + r1) * s[1]}};
    }

   private:
    [[nodiscard]] stl::vector<real_t> find_peak_positions_for_plane(const xrd::crystal& c, const rvec3_t& plane) const {
      const xrd::single_plane_diffraction_pattern experiment(c, m_Size, m_MosaicSpread, 1000, plane, m_Temperature, m_Wavelength);

      const auto intensities = experiment.generate(m_Angles);
      auto peak_indices = math::find_peak_indices(intensities);
      std::sort(peak_indices.begin(), peak_indices.end());

      stl::vector<real_t> peaks(peak_indices.size());
      std::transform(peak_indices.begin(), peak_indices.end(), peaks.begin(), [in = std::span(m_Angles)](sint_t i) { return in[i]; });

      return peaks;
    }

    [[nodiscard]] real_t secondary_peaks_energy(const xrd::crystal& c) const {
      //      auto peaks_110 = find_peak_positions_for_plane(c, {1, 1, 0});
      //      if(peaks_110.size() < 1)
      //        return std::numeric_limits<real_t>::max();
      //      real_t e_110 = math::sqr(m_Peak_111 - 2 * peaks_110[0]);
      //
      //      auto peaks_111 = find_peak_positions_for_plane(c, {1, 1, 1});
      //      if(peaks_111.size() < 1)
      //        return std::numeric_limits<real_t>::max();
      //      real_t e_111 = math::sqr(m_Peak_200 - 2 * peaks_111[0]);
      //
      //      auto peaks_200 = find_peak_positions_for_plane(c, {2, 0, 0});
      //      if(peaks_200.size() < 1)
      //        return std::numeric_limits<real_t>::max();
      //      real_t e_100 = math::sqr(m_Peak_110 - 2 * peaks_200[0]);

      auto peaks_110 = find_peak_positions_for_plane(c, {1, 1, 0});
      if(peaks_110.empty())
        return std::numeric_limits<real_t>::max();
      real_t e_110 = math::sqr(m_Peak_110 - 2 * peaks_110[0]);

      auto peaks_111 = find_peak_positions_for_plane(c, {1, 1, 1});
      if(peaks_111.empty())
        return std::numeric_limits<real_t>::max();
      real_t e_111 = math::sqr(m_Peak_111 - 2 * peaks_111[0]);

      auto peaks_200 = find_peak_positions_for_plane(c, {2, 0, 0});
      if(peaks_200.size() < 2)
        return std::numeric_limits<real_t>::max();
      real_t e_100 = math::sqr(m_Peak_200 - 2 * peaks_200[1]);

      //      fmt::print("[{}]\n", fmt::join(peaks_200, ", "));

      return e_110 + e_111 + e_100;
    }

    ivector_t<3> m_Size;
    real_t m_MosaicSpread;

    real_t m_Temperature;
    real_t m_Wavelength;

    rdata_t m_Angles;

    real_t m_Peak_001, m_Peak_110, m_Peak_111, m_Peak_200, m_Peak_002;
  };
}    // namespace

template <>
struct fmt::formatter<xrd_annealing_simulation::solution_type> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(const xrd_annealing_simulation::solution_type& s, FormatContext& ctx) {
    return format_to(ctx.out(), "(a = {:.5f}; c = {:.5f})", s[0], s[1]);
  }
};

int main(int argc, char** argv) {
  fmt::format("{}", xrd_annealing_simulation::solution_type{});

  hr_timer timer{"Annealing"};

  opt::simulated_annealer<xrd_annealing_simulation, true> annealer;
  const xrd_annealing_simulation xas({40, 40, 40}, math::deg2rad(0.5), 300, xray::CuKalpha::lambda, math::data::linspace(10, 30, 2000));

  timer.start();
  auto s = annealer.run(xas, 10, 100);
  timer.stop();
  timer.report();

  fmt::print("[{}]\n", fmt::join(s, ", "));
}

//int main(int argc, char** argv) {
//  ds::dataset_2d pattern(io::load_csv("fept/AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta.txt"));
//  pattern.y() = pattern.y().log();
//
//  auto [peak_pos, peak_mag] = pattern.find_peaks_in_interval(20, 55, 0.1);
//
//  fmt::print("[{}]\n", fmt::join(peak_pos, ", "));
//}