#include "convolution.hpp"

math::convolution_kernel_1d math::convolution_kernel_1d::boxcar(sint_t n) noexcept {
  convolution_kernel_1d ck(n);
  ck.modify([&ck](const sint_t ii) -> real_t { return 1 / real_t(ck.size()); });
  return ck;
}

math::convolution_kernel_1d math::convolution_kernel_1d::gaussian(real_t sigma) noexcept {
  convolution_kernel_1d ck(std::ceil(3 * sigma));
  ck.modify([c1 = 1 / std::sqrt(2 * C_PI * sigma * sigma), c2 = 1 / (2 * sigma * sigma)](const sint_t ii) -> real_t { return c1 * std::exp(-c2 * ii * ii); });
  return ck;
}

rdata_t math::convolution_kernel_1d::apply(std::span<const real_t> signal) const {
  rdata_t buf(signal.size());
  {
    auto blurred = std::span(buf);
    auto fn_get_partial_block = [this, signal](const sint_t i) -> std::span<const real_t> {
      const sint_t l_bound = std::max<sint_t>(0, i - m_N), u_bound = std::min<sint_t>(signal.size(), i + m_N + 1);
      return signal.subspan(l_bound, u_bound - l_bound);
    };
    auto fn_get_full_block = [this, signal](const sint_t i) -> std::span<const real_t> {
      return signal.subspan(i - m_N, 2 * m_N + 1);
    };

    if(static_cast<sint_t>(blurred.size()) >= (2 * m_N + 1)) {
      auto l_bound = m_N;
      auto r_bound = static_cast<sint_t>(blurred.size()) - m_N - 1;

      for(sint_t ii = 0; ii < l_bound; ++ii)
        blurred[ii] = apply_partial(fn_get_partial_block(ii), ii);
      for(sint_t ii = l_bound; ii < r_bound; ++ii)
        blurred[ii] = apply_full(fn_get_full_block(ii));
      for(sint_t ii = r_bound; ii < static_cast<sint_t>(blurred.size()); ++ii)
        blurred[ii] = apply_partial(fn_get_partial_block(ii), ii - r_bound);
    } else {
      for(sint_t ii = 0; ii < static_cast<sint_t>(blurred.size()); ++ii)
        blurred[ii] = apply_partial(fn_get_partial_block(ii), ii);
    }
  }
  return buf;
}
