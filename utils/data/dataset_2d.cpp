#include "dataset_2d.hpp"

#include <concepts/compare.hpp>
#include <numeric>
#include <span>

#include <fmt/format.h>

#include "math/peak_finder.hpp"

namespace {
  inline ds::dataset_2d find_peaks(std::span<const real_t> x, std::span<const real_t> y, real_t relative_selectivity, real_t threshold,
                                   math::peak_type extrema) {
    auto peak_ind = math::find_peak_indices(y, relative_selectivity, threshold, extrema);

    stl::vector<real_t> peak_pos, peak_mag;
    peak_pos.resize(peak_ind.size());
    std::transform(peak_ind.begin(), peak_ind.end(), peak_pos.begin(), [x](sint_t i) { return x[i]; });
    peak_mag.resize(peak_ind.size());
    std::transform(peak_ind.begin(), peak_ind.end(), peak_mag.begin(), [y](sint_t i) { return y[i]; });

    return {std::move(peak_pos), std::move(peak_mag)};
  }

  template <std::totally_ordered T>
  inline std::tuple<std::span<T>, std::span<T>> get_interval(std::span<T> x, std::span<T> y, std::remove_const_t<T> x_min, std::remove_const_t<T> x_max) {
    if(x_min == x_max)
      throw std::invalid_argument(fmt::format("invalid interval: [{}; {}]", x_min, x_max));
    if(x_min > x_max)
      std::swap(x_min, x_max);

    auto x_min_ind = std::distance(x.begin(), std::lower_bound(x.begin(), x.end(), x_min));
    auto x_max_ind = std::distance(x.begin(), std::upper_bound(x.begin(), x.end(), x_max));

    auto span_len = x_max_ind - x_min_ind;
    if(span_len == 0)
      throw std::invalid_argument(fmt::format("invalid interval: [{}; {}] is empty", x_min, x_max));

    return {x.subspan(x_min_ind, span_len), y.subspan(x_min_ind, span_len)};
  }

  inline const real_t& get_value(std::span<const real_t> x, std::span<const real_t> y, real_t x_val) {
    auto val_it = std::lower_bound(x.begin(), x.end(), x_val);
    if(val_it == x.end())
      throw std::runtime_error(fmt::format("out of range: {0} not in [{1}, {2}]", x_val, x.front(), x.back()));
    if(*val_it != x_val)
      throw std::runtime_error(fmt::format("out of range: {0} not an element in x set", x_val));

    return y[std::distance(x.begin(), val_it)];
  }

  inline real_t get_interpolated_value(std::span<const real_t> x, std::span<const real_t> y, real_t x_val) {
    auto val_it = std::lower_bound(x.begin(), x.end(), x_val);
    if(val_it == x.end())
      throw std::runtime_error(fmt::format("out of range: {0} not in [{1}, {2}]", x_val, x.front(), x.back()));
    if(*val_it == x_val)
      return y[std::distance(x.begin(), val_it)];
    if(val_it == x.begin())
      throw std::runtime_error(fmt::format("out of range: {0} not in [{1}, {2}]", x_val, x.front(), x.back()));

    auto val_ind = std::distance(x.begin(), val_it);
    auto x_1 = x[val_ind - 1], y_1 = y[val_ind - 1];
    auto m = (y[val_ind] - y_1) / (x[val_ind] - x_1);

    return m * (x_val - x_1) + y_1;
  }

  inline void validate_dataset(std::span<const real_t> x, std::span<const real_t> y) {
    if(x.size() != y.size())
      throw std::invalid_argument(fmt::format("shape mismatch: x and y have different sizes ({} != {})", x.size(), y.size()));
    if(!std::is_sorted(x.begin(), x.end()))
      throw std::invalid_argument("x is not sorted");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
      throw std::invalid_argument("x contains duplicate values");
  }
}    // namespace

ds::dataset_2d::dataset_2d(stl::vector<real_t> x, stl::vector<real_t> y) {
  if(x.size() != y.size())
    throw std::invalid_argument(fmt::format("shape mismatch: x and y have different sizes ({} != {})", x.size(), y.size()));

  if(!std::is_sorted(x.begin(), x.end())) {
    stl::vector<uint_t> sorted_indices(x.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), uint_t(0));
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&x](uint_t i1, uint_t i2) { return x[i1] < x[i2]; });

    for(uint_t ii = 0; ii < sorted_indices.size(); ++ii) {
      std::swap(x[ii], x[sorted_indices[ii]]);
      std::swap(y[ii], y[sorted_indices[ii]]);
      std::swap(sorted_indices[ii], sorted_indices[sorted_indices[ii]]);
    }
  }

  ::validate_dataset(x, y);

  m_X = std::move(x);
  m_Y = std::move(y);
}

ds::dataset_2d_span::dataset_2d_span(std::span<real_t> x, std::span<real_t> y) : m_X(x), m_Y(y) {
  ::validate_dataset(m_X, m_Y);
}

ds::dataset_2d ds::dataset_2d_span::find_peaks(real_t relative_selectivity, real_t threshold, math::peak_type extrema) const {
  if(m_X.size() < 3)
    throw std::runtime_error("dataset is too small");

  return ::find_peaks(m_X, m_Y, relative_selectivity, threshold, extrema);
}

ds::dataset_2d_span ds::dataset_2d_span::get(real_t x_min, real_t x_max) const {
  auto [x_span, y_span] = get_interval(m_X, m_Y, x_min, x_max);
  return {x_span, y_span, no_validation};
}

ds::dataset_2d_view::dataset_2d_view(std::span<const real_t> x, std::span<const real_t> y) : m_X(x), m_Y(y) {
  ::validate_dataset(m_X, m_Y);
}

ds::dataset_2d ds::dataset_2d_view::find_peaks(real_t relative_selectivity, real_t threshold, math::peak_type extrema) const {
  return ::find_peaks(m_X, m_Y, relative_selectivity, threshold, extrema);
}

real_t ds::dataset_2d_view::get(real_t x) const {
  return get_interpolated_value(m_X, m_Y, x);
}

real_t ds::dataset_2d_view::get(real_t x, no_interpolation_tag) const {
  return get_value(m_X, m_Y, x);
}

ds::dataset_2d_view ds::dataset_2d_view::get(real_t x_min, real_t x_max) const {
  auto [x_span, y_span] = get_interval(m_X, m_Y, x_min, x_max);
  return {x_span, y_span, no_validation};
}
