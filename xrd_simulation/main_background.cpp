#include <data/dataset_2d.hpp>
#include <io.hpp>
#include <math/convolution.hpp>
#include <optimisation/simulated_annealing.hpp>
#include <timer.hpp>

namespace {
  class background_processing {
    struct interval {
      real_t x_min, x_max;
    };

   public:
    background_processing(std::string_view measurement_path, std::string_view background_path)
        : m_Measurement{io::load_csv(measurement_path)}, m_Background{io::load_csv(background_path)} {}

    std::tuple<ds::dataset_2d, ds::dataset_2d> process() const;

   private:
    ds::dataset_2d m_Measurement, m_Background;
  };
}    // namespace

std::tuple<ds::dataset_2d, ds::dataset_2d> background_processing::process() const {
  constexpr std::array<interval, 3> fitting_intervals = {{{27, 31}, {61, 67}, {105, 110}}};

  std::array<real_t, fitting_intervals.size()> bg_means, m_means;
  std::transform(fitting_intervals.begin(), fitting_intervals.end(), bg_means.begin(),
                 [this](const interval& i) { return m_Background.get(i.x_min, i.x_max).arithmetic_mean(); });
  std::transform(fitting_intervals.begin(), fitting_intervals.end(), m_means.begin(),
                 [this](const interval& i) { return m_Measurement.get(i.x_min, i.x_max).arithmetic_mean(); });

  std::array<real_t, fitting_intervals.size()> interval_factors;
  std::transform(bg_means.begin(), bg_means.end(), m_means.begin(), interval_factors.begin(), [](real_t bg, real_t m) { return m / bg; });

  real_t factor = std::reduce(interval_factors.begin(), interval_factors.end()) / fitting_intervals.size();

  ds::dataset_2d new_bg = m_Background, new_m = m_Measurement;
  new_bg.y() = factor * math::convolution_kernel_1d::boxcar(5).apply(m_Background.y());
  new_m.y() -= new_bg.y();

  fmt::print("Correction factor: {0}\n", factor);
  return std::make_tuple<ds::dataset_2d, ds::dataset_2d>(std::move(new_m), std::move(new_bg));
}

int main() {
  auto bg_processing = background_processing("fept/AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta.txt", "fept/MgO_Background_XRD_20_120_2-Theta_Omega.txt");

  ds::dataset_2d processed_measurement, processed_background;
  benchmark("Background processing", [&bg_processing, &processed_measurement, &processed_background]() {
    std::tie(processed_measurement, processed_background) = bg_processing.process();
  });

  io::write_csv("fept/AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta_signal.txt", processed_measurement);
  io::write_csv("fept/MgO_Background_XRD_20_120_2-Theta_Omega_processed.txt", processed_background);
}
