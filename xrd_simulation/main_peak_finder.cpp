#include <fmt/format.h>

#include <data/dataset_2d.hpp>
#include <io.hpp>
#include <math.hpp>
#include <timer.hpp>
#include <types_format.hpp>

template <>
struct fmt::formatter<ds::dataset_2d_view::point> {
  auto parse(format_parse_context& ctx) {
    return ctx.begin();
  }

  auto format(const ds::dataset_2d_view::point& p, format_context& ctx) {
    return format_to(ctx.out(), "({0}, {1})", p.x, p.y);
  }
};

int main() {
//  real_t angle_min = 52.5, angle_max = 60;
//
//  ds::dataset_2d simulated_signal(io::load_csv("../../Report/S1/fig/fept_measurement/xrd_full_normalised.csv"));
//  ds::dataset_2d measured_signal(io::load_csv("../../Report/S1/fig/fept_measurement/AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta_signal.txt"));
//
//  fmt::print("{{{}}}\n", fmt::join(simulated_signal.get(angle_min, angle_max).find_peaks(), ", "));
//  fmt::print("{{{}}}\n", fmt::join(measured_signal.get(2*angle_min, 2*angle_max).find_peaks(0.5), ", "));

  ds::dataset_2d signal_1(io::load_csv("fept/MgO_Background_XRD_20_120_2-Theta_Omega.txt"));
  ds::dataset_2d signal_2(io::load_csv("fept/AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta.txt"));

  fmt::print("{{{}}}\n", fmt::join(signal_1.find_peaks(), ", "));
  fmt::print("{{{}}}\n", fmt::join(signal_2.find_peaks(0.5), ", "));
}