#ifndef XRD_PEAK_FINDER_HPP
#define XRD_PEAK_FINDER_HPP

#include <vector>

#include <span>

#include <types.hpp>

namespace math {
  //  stl::vector<sint_t> find_peaks(std::span<const real_t> x0, real_t relative_selectivity = 0.25);

  enum class peak_type { e_Minima = -1, e_Maxima = 1 };

  stl::vector<sint_t> find_peak_indices(std::span<const real_t> x0, real_t relative_selectivity = 0.25, real_t threshold = 0,
                                        peak_type extrema = peak_type::e_Maxima, bool include_endpoints = false);

  inline std::pair<stl::vector<sint_t>, stl::vector<real_t>> find_peaks(std::span<const real_t> x0, real_t relative_selectivity = 0.25, real_t threshold = 0,
                                                                        peak_type extrema = peak_type::e_Maxima, bool include_endpoints = false) {
    std::pair<stl::vector<sint_t>, stl::vector<real_t>> peaks;

    peaks.first = find_peak_indices(x0, relative_selectivity, threshold, extrema, include_endpoints);

    peaks.second.resize(peaks.first.size());
    std::transform(peaks.first.begin(), peaks.first.end(), peaks.second.begin(), [x0](sint_t i) { return x0[i]; });

    return peaks;
  }
}    // namespace math

#endif    //XRD_PEAK_FINDER_HPP
