#ifndef XRD_DATASET_2D_HPP
#define XRD_DATASET_2D_HPP

#include <range/v3/iterator/basic_iterator.hpp>

#include "math/peak_finder.hpp"
#include "types.hpp"

namespace ds {
  namespace detail {
    template <typename ElementType, typename ValueType, typename ReadType = ValueType>
    struct iterator_base {
      using span_iterator_type = typename std::span<ElementType>::iterator;

      using value_type = ValueType;

      struct mixin;

      iterator_base() = default;
      iterator_base(span_iterator_type x_iter, span_iterator_type y_iter) : x_iterator{x_iter}, y_iterator{y_iter} {}

      [[nodiscard]] ReadType read() const {
        return {*x_iterator, *y_iterator};
      }

      [[nodiscard]] bool equal(const iterator_base& other) const {
        return x_iterator == other.x_iterator;
      }

      void next() {
        ++x_iterator;
        ++y_iterator;
      }
      void prev() {
        --x_iterator;
        --y_iterator;
      }
      void advance(std::ptrdiff_t n) {
        x_iterator += n;
        y_iterator += n;
      }

      [[nodiscard]] std::ptrdiff_t distance_to(const iterator_base& other) const {
        return other.x_iterator - this->x_iterator;
      }

     private:
      span_iterator_type x_iterator, y_iterator;
    };

    template <typename ElementType, typename ValueType, typename ReadType>
    struct iterator_base<ElementType, ValueType, ReadType>::mixin : ranges::basic_mixin<iterator_base<ElementType, ValueType, ReadType>> {
      using ranges::basic_mixin<iterator_base<ElementType, ValueType, ReadType>>::basic_mixin;

      inline mixin(span_iterator_type x_iter, span_iterator_type y_iter) : mixin{iterator_base{x_iter, y_iter}} {}
    };

    template <typename ElementType, typename ValueType, typename ReadType = ValueType>
    using iterator = ranges::basic_iterator<iterator_base<ElementType, ValueType, ReadType>>;
  }

  struct no_interpolation_tag {};
  inline constexpr no_interpolation_tag no_interpolation;

  struct no_validation_tag {};
  inline constexpr no_validation_tag no_validation;

  class dataset_2d;

  class dataset_2d_view {
   public:
    struct point {
      real_t x, y;
    };

   private:
    using iterator = detail::iterator<const real_t, point>;

   public:
    inline dataset_2d_view() noexcept = default;
    dataset_2d_view(std::span<const real_t> x, std::span<const real_t> y);
    inline dataset_2d_view(std::span<const real_t> x, std::span<const real_t> y, no_validation_tag) noexcept : m_X(x), m_Y(y) {}
    inline explicit dataset_2d_view(std::tuple<std::span<const real_t>, std::span<const real_t>> dset)
        : dataset_2d_view(std::get<0>(dset), std::get<1>(dset)) {}

    [[nodiscard]] inline auto begin() const noexcept {
      return iterator{m_X.begin(), m_Y.begin()};
    }
    [[nodiscard]] inline auto end() const noexcept {
      return iterator{m_X.end(), m_Y.end()};
    }

    [[nodiscard]] inline real_t arithmetic_mean() const noexcept {
      return y().mean();
    }
    [[nodiscard]] dataset_2d find_peaks(real_t relative_selectivity = 0.25, real_t threshold = 0, math::peak_type extrema = math::peak_type::e_Maxima) const;

    [[nodiscard]] real_t get(real_t x) const;
    [[nodiscard]] real_t get(real_t x, no_interpolation_tag) const;
    [[nodiscard]] dataset_2d_view get(real_t x_min, real_t x_max) const;

    [[nodiscard]] inline rdata_view_t x() const noexcept {
      return rdata_view_t(m_X.data(), m_X.size());
    }
    [[nodiscard]] inline rdata_view_t y() const noexcept {
      return rdata_view_t(m_Y.data(), m_Y.size());
    }

    [[nodiscard]] inline size_t size() const noexcept {
      return m_X.size();
    }

    [[nodiscard]] inline point operator[](size_t index) const noexcept {
      return {m_X[index], m_Y[index]};
    }

   private:
    std::span<const real_t> m_X, m_Y;
  };

  class dataset_2d_span {
    using iterator = detail::iterator<real_t, dataset_2d_view::point>;

   public:
    inline dataset_2d_span() noexcept = default;
    dataset_2d_span(std::span<real_t> x, std::span<real_t> y);
    inline dataset_2d_span(std::span<real_t> x, std::span<real_t> y, no_validation_tag) noexcept : m_X(x), m_Y(y) {}
    inline explicit dataset_2d_span(std::tuple<std::span<real_t>, std::span<real_t>> dset) : dataset_2d_span(std::get<0>(dset), std::get<1>(dset)) {}

    [[nodiscard]] inline auto begin() const noexcept {
      return iterator{m_X.begin(), m_Y.begin()};
    }
    [[nodiscard]] inline auto end() const noexcept {
      return iterator{m_X.end(), m_Y.end()};
    }

    [[nodiscard]] inline real_t arithmetic_mean() const noexcept {
      return dataset_2d_view(*this).arithmetic_mean();
    }
    [[nodiscard]] dataset_2d find_peaks(real_t relative_selectivity = 0.25, real_t threshold = 0, math::peak_type extrema = math::peak_type::e_Maxima) const;

    [[nodiscard]] inline real_t get(real_t x) const {
      return dataset_2d_view(*this).get(x);
    }
    [[nodiscard]] inline real_t get(real_t x, no_interpolation_tag) const {
      return dataset_2d_view(*this).get(x, no_interpolation);
    }
    [[nodiscard]] dataset_2d_span get(real_t x_min, real_t x_max) const;

    [[nodiscard]] inline rdata_view_t x() const noexcept {
      return rdata_view_t(m_X.data(), m_X.size());
    }
    [[nodiscard]] inline rdata_span_t y() const noexcept {
      return rdata_span_t(m_Y.data(), m_Y.size());
    }

    [[nodiscard]] inline size_t size() const noexcept {
      return m_X.size();
    }

    [[nodiscard]] inline auto operator[](size_t index) const noexcept {
      return dataset_2d_view(*this)[index];
    }

    inline operator dataset_2d_view() const noexcept {
      return {m_X, m_Y, no_validation};
    }

   private:
    std::span<real_t> m_X, m_Y;
  };

  class dataset_2d {
    using dataset_type = std::pair<stl::vector<real_t>, stl::vector<real_t>>;

   public:
    inline dataset_2d() noexcept = default;
    inline explicit dataset_2d(dataset_2d_view dset)
        : dataset_2d(stl::vector<real_t>(dset.x().data(), dset.x().data() + dset.x().size()),
                     stl::vector<real_t>(dset.y().data(), dset.y().data() + dset.y().size()), no_validation) {}
    dataset_2d(stl::vector<real_t> x, stl::vector<real_t> y);
    inline explicit dataset_2d(std::tuple<stl::vector<real_t>, stl::vector<real_t>> dset)
        : dataset_2d(std::move(std::get<0>(dset)), std::move(std::get<1>(dset))) {}
    inline dataset_2d(stl::vector<real_t> x, stl::vector<real_t> y, no_validation_tag) : m_X(std::move(x)), m_Y(std::move(y)) {}
    inline dataset_2d(const dataset_2d&) = default;
    inline dataset_2d(dataset_2d&&) noexcept = default;

    inline dataset_2d& operator=(const dataset_2d&) = default;
    inline dataset_2d& operator=(dataset_2d&&) noexcept = default;

    [[nodiscard]] inline auto begin() noexcept {
      return dataset_2d_span(*this).begin();
    }
    [[nodiscard]] inline auto begin() const noexcept {
      return dataset_2d_view(*this).begin();
    }
    [[nodiscard]] inline auto end() noexcept {
      return dataset_2d_span(*this).end();
    }
    [[nodiscard]] inline auto end() const noexcept {
      return dataset_2d_view(*this).end();
    }

    [[nodiscard]] inline real_t arithmetic_mean() const noexcept {
      return dataset_2d_view(*this).arithmetic_mean();
    }
    [[nodiscard]] inline dataset_2d find_peaks(real_t relative_selectivity = 0.25, real_t threshold = 0,
                                               math::peak_type extrema = math::peak_type::e_Maxima) const {
      return dataset_2d_view(*this).find_peaks(relative_selectivity, threshold, extrema);
    }

    [[nodiscard]] inline real_t get(real_t x) const {
      return dataset_2d_view(*this).get(x);
    }
    [[nodiscard]] inline real_t get(real_t x, no_interpolation_tag) const {
      return dataset_2d_view(*this).get(x, no_interpolation);
    }
    [[nodiscard]] inline dataset_2d_span get(real_t x_min, real_t x_max) {
      return dataset_2d_span(*this).get(x_min, x_max);
    }
    [[nodiscard]] inline dataset_2d_view get(real_t x_min, real_t x_max) const {
      return dataset_2d_view(*this).get(x_min, x_max);
    }

    [[nodiscard]] inline rdata_view_t x() const noexcept {
      return rdata_view_t(m_X.data(), m_X.size());
    }
    [[nodiscard]] inline rdata_span_t y() noexcept {
      return rdata_span_t(m_Y.data(), m_Y.size());
    }
    [[nodiscard]] inline rdata_view_t y() const noexcept {
      return rdata_view_t(m_Y.data(), m_Y.size());
    }

    [[nodiscard]] inline size_t size() const noexcept {
      return m_X.size();
    }

    [[nodiscard]] inline auto operator[](size_t index) const noexcept {
      return dataset_2d_view(*this)[index];
    }

    inline operator dataset_2d_span() noexcept {
      return {m_X, m_Y, no_validation};
    }
    inline operator dataset_2d_view() const noexcept {
      return {m_X, m_Y, no_validation};
    }

   private:
    stl::vector<real_t> m_X, m_Y;
  };
}    // namespace ds

#endif    //XRD_DATASET_2D_HPP
