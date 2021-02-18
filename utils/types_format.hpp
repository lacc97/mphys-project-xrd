#ifndef XRD_TYPES_FORMAT_HPP
#define XRD_TYPES_FORMAT_HPP

#include <fmt/format.h>

#include "types.hpp"

template <typename Numeric>
struct fmt::formatter<data_t<Numeric>> {
  auto parse(format_parse_context& ctx) {
    return ctx.begin();
  }

  auto format(const data_t<Numeric>& d, format_context& ctx) {
    return fmt::format_to(ctx.out(), "[{}]", fmt::join(d.data(), d.data() + d.size(), ", "));
  }
};

template <typename Numeric>
struct fmt::formatter<data_span_t<Numeric>> {
  auto parse(format_parse_context& ctx) {
    return ctx.begin();
  }

  auto format(const data_span_t<Numeric>& d, format_context& ctx) {
    return fmt::format_to(ctx.out(), "[{}]", fmt::join(d.data(), d.data() + d.size(), ", "));
  }
};

template <typename Numeric>
struct fmt::formatter<data_view_t<Numeric>> {
  auto parse(format_parse_context& ctx) {
    return ctx.begin();
  }

  auto format(const data_view_t<Numeric>& d, format_context& ctx) {
    return fmt::format_to(ctx.out(), "[{}]", fmt::join(d.data(), d.data() + d.size(), ", "));
  }
};

template <typename Numeric, int N>
struct fmt::formatter<vector_t<Numeric, N>> {
  auto parse(format_parse_context& ctx) {
    return ctx.begin();
  }

  auto format(const vector_t<Numeric, N>& v, format_context& ctx) {
    return fmt::format_to(ctx.out(), "({})", fmt::join(v.data(), v.data() + v.size(), ", "));
  }
};

#endif    //XRD_TYPES_FORMAT_HPP
