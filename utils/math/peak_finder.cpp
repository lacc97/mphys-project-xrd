#include "peak_finder.hpp"

#include <array>

#include <allocator.hpp>
#include <math.hpp>

/// Based on: https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate
/// Based on: https://github.com/claydergc/find-peaks

namespace {
  inline void diff(std::span<const real_t> in, std::span<real_t> out) {
    auto len = std::min(in.size(), out.size() + 1);
    for(uint_t ii = 1; ii < len; ++ii)
      out[ii - 1] = in[ii] - in[ii - 1];
  }
  inline stl::vector<real_t> diff(std::span<const real_t> in) {
    stl::vector<real_t> out(in.size() - 1);
    diff(in, out);
    return out;
  }

  inline stl::vector<sint_t> find_indices_less_than(std::span<const real_t> in, real_t threshold) {
    stl::vector<sint_t> out;
    for(uint_t i = 0; i < in.size(); ++i)
      if(in[i] < threshold)
        out.push_back(i + 1);
    return out;
  }
}    // namespace

stl::vector<sint_t> math::find_peak_indices(std::span<const real_t> x0, real_t relative_selectivity, real_t threshold, math::peak_type extrema,
                                            bool include_endpoints) {
  stl::vector<real_t> x0_buf;
  if(extrema == peak_type::e_Minima) {
    /* make it as if we are finding maxima regardless */
    x0_buf.resize(x0.size());
    std::transform(x0.begin(), x0.end(), x0_buf.begin(), [](real_t r) -> real_t { return -r; });
    x0 = x0_buf;
  }

  /* calculate absolute selectivity */
  real_t selectivity;
  {
    auto min_index = std::distance(x0.begin(), std::min_element(x0.begin(), x0.end()));
    auto max_index = std::distance(x0.begin(), std::max_element(x0.begin(), x0.end()));

    selectivity = (x0[max_index] - x0[min_index]) * relative_selectivity;
  }

  /* adjust threshold according to extrema */
  threshold *= (extrema == peak_type::e_Maxima ? 1 : -1);

  stl::vector<sint_t> ind;
  {
    /* find derivative */
    stl::vector<real_t> dx = diff(x0);
    std::replace(dx.begin(), dx.end(), real_t(0.0), -std::numeric_limits<real_t>::epsilon());

    std::span<real_t> dx0 = std::span(dx).first(dx.size() - 1);
    std::span<real_t> dx1 = std::span(dx).last(dx.size() - 1);

    stl::vector<real_t> dx2(dx.size() - 1);
    std::transform(dx0.begin(), dx0.end(), dx1.begin(), dx2.begin(), std::multiplies<>());

    /* find where the derivative changes sign */
    ind = find_indices_less_than(dx2, 0);
  }

  stl::vector<real_t> x;
  real_t min_mag, left_min;
  if(include_endpoints) {
    x.resize(ind.size() + 2);
    x[0] = x0[0];
    std::transform(ind.begin(), ind.end(), x.begin() + 1, [x0](sint_t i) -> real_t { return x0[i]; });
    x[x.size() - 1] = x0[x0.size() - 1];

    ind.insert(ind.begin(), 0);
    ind.insert(ind.end(), x0.size());

    min_mag = *std::min_element(x.begin(), x.end());
    left_min = min_mag;
  } else {
    x.resize(ind.size());
    std::transform(ind.begin(), ind.end(), x.begin(), [x0](sint_t i) -> real_t { return x0[i]; });

    min_mag = *std::min_element(x.begin(), x.end());
    left_min = std::min(x[0], x0[0]);
  }

  stl::vector<sint_t> peak_indices;
  if(x.size() > 2) {
    real_t temp_mag = min_mag;
    bool found_peak = false;
    uint_t ii;

    if(include_endpoints) {
      /* Deal with first point a little differently since we tacked it on.
     * Calculate the sign of the derivative since we tacked the first point on
     *  it does not neccessarily alternate like the rest.*/

      std::array<real_t, 2> sign_dx;
      {
        std::array<real_t, 2> x_diff;
        diff(std::span(x).first(3), x_diff);

        std::transform(x_diff.begin(), x_diff.end(), sign_dx.begin(), math::signum<real_t>);
      }

      /* we want alternating signs */
      if(sign_dx[0] == sign_dx[1]) {
        if(sign_dx[0] <= 0) {
          x.erase(x.begin() + 1);
          ind.erase(ind.begin() + 1);
        } else {
          x.erase(x.begin());
          ind.erase(ind.begin());
        }
      }
    }

    /* skip the first point if it is smaller so we always start on a maxima */
    if(x[0] >= x[1])
      ii = 0;
    else
      ii = 1;

    /* preallocate max number of maxima */
    stl::vector<sint_t> peak_loc;
    stl::vector<real_t> peak_mag;
    {
      uint_t max_peaks = std::ceil(x.size() / 2.0);

      peak_loc = stl::vector<sint_t>(max_peaks, 0);
      peak_mag = stl::vector<real_t>(max_peaks, 0.0);
    }

    sint_t c_index = 1;
    sint_t temp_loc;

    /* loop through extrema, which should be peaks and then valleys */
    while(ii < x.size()) {
      /* start with a peak */
      ii = ii + 1;

      /* reset peak finding if we already had a peak and the next peak is bigger
       *  than the last or the left min was small enough to reset */
      if(found_peak) {
        temp_mag = min_mag;
        found_peak = false;
      }

      /* found new peak that was larger than temp_mag and selectivity larger
       *  than the minimum to its left.*/
      if(x[ii - 1] > temp_mag && x[ii - 1] > left_min + selectivity) {
        temp_loc = ii - 1;
        temp_mag = x[ii - 1];
      }

      /* make sure we don't iterate past the length of our vector
       * we assign the last point differently out of the loop */
      if(ii == x.size())
        break;

      /* move onto valley */
      ii = ii + 1;

      /* come down at least sel from peak */
      if(temp_mag > selectivity + x[ii - 1]) {
        /* we have found a peak */
        found_peak = true;
        left_min = x[ii - 1];
        peak_loc[c_index - 1] = temp_loc;
        peak_mag[c_index - 1] = temp_mag;
        c_index = c_index + 1;
      } else if(x[ii - 1] < left_min) {
        /* we have a new left minima */
        left_min = x[ii - 1];
      }
    }

    /* check end point */
    if(x[x.size() - 1] > temp_mag && x[x.size() - 1] > left_min + selectivity) {
      found_peak = true;
      peak_loc[c_index - 1] = x.size() - 1;
      peak_mag[c_index - 1] = x[x.size() - 1];
      c_index = c_index + 1;
    }
    if(!found_peak) {
      /* check if we still have to add the last point */
      if(include_endpoints) {
        if(temp_mag > min_mag) {
          peak_loc[c_index - 1] = temp_loc;
          peak_mag[c_index - 1] = temp_mag;
          c_index = c_index + 1;
        }
      } else {
        if(temp_mag > std::min(x.back(), x0.back()) + selectivity) {
          peak_loc[c_index - 1] = temp_loc;
          peak_mag[c_index - 1] = temp_mag;
          c_index = c_index + 1;
        }
      }
    }

    /* Create output */
    if(c_index > 0) {
      std::span<sint_t> peak_loc_temp = std::span(peak_loc).first(c_index - 1);
      std::transform(peak_loc_temp.begin(), peak_loc_temp.end(), std::back_inserter(peak_indices), [in = std::span(ind)](sint_t i) { return in[i]; });
    }
    return peak_indices;
  } else {
    /* this is a monotone function where the endpoint is the only peak */
    auto x_max = std::max_element(x.begin(), x.end());
    if(include_endpoints && *x_max > min_mag + selectivity)
      peak_indices.push_back(std::distance(x.begin(), x_max));
  }

  /* apply threshold */
  peak_indices.erase(std::remove_if(peak_indices.begin(), peak_indices.end(), [threshold](real_t m) -> bool { return m > threshold; }), peak_indices.end());

  return peak_indices;
}
