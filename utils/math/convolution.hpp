#ifndef XRD_CONVOLUTION_HPP
#define XRD_CONVOLUTION_HPP

#include <numeric>
#include <span>

#include "constants.hpp"
#include "types.hpp"

namespace math {
  class convolution_kernel_1d {
   public:
    static convolution_kernel_1d boxcar(sint_t n) noexcept;
    static convolution_kernel_1d gaussian(real_t ii) noexcept;

    rdata_t apply(std::span<const real_t> signal) const;

    [[nodiscard]] inline sint_t n() const noexcept {
      return m_N;
    }
    [[nodiscard]] inline sint_t size() const noexcept {
      return m_Kernel.size();
    }

    inline real_t operator[](sint_t index) const noexcept {
      return m_Kernel(index + m_N);
    }

   private:
    convolution_kernel_1d(sint_t n) : m_N{n}, m_Kernel{2 * m_N + 1} {}

    inline real_t apply_full(std::span<const real_t> block) const noexcept {
      return m_Kernel.dot(rvector_view_t<n_dynamic>(block.data(), block.size()));
    }
    inline real_t apply_partial(std::span<const real_t> block, sint_t index) const noexcept {
      std::span<const real_t> kernel_span;
      {
        const sint_t l_bound = std::max<sint_t>(0, index - static_cast<sint_t>(block.size())), u_bound = std::min<sint_t>(block.size(), index + m_N + 1);
        kernel_span = std::span(m_Kernel).subspan(l_bound, u_bound - l_bound);
      }

      auto kern = rvector_view_t<n_dynamic>(kernel_span.data(), kernel_span.size());
      return kern.dot(rvector_view_t<n_dynamic>(block.data(), block.size())) / kern.sum();
    }

    template <typename F>
    inline void modify(F&& f) noexcept {
      for(sint_t ii = -m_N; ii <= m_N; ++ii)
        (*this)[ii] = f(ii);
      m_Kernel /= m_Kernel.sum();
    }

    inline real_t& operator[](sint_t index) noexcept {
      return m_Kernel(index + m_N);
    }

    sint_t m_N;
    rvector_t<n_dynamic> m_Kernel;
  };
}    // namespace math

#endif    //XRD_CONVOLUTION_HPP
