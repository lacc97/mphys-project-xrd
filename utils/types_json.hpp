#ifndef XRD_TYPES_JSON_HPP
#define XRD_TYPES_JSON_HPP

#include <nlohmann/json.hpp>

namespace nlohmann {
  template <typename Numeric>
  struct adl_serializer<data_view_t<Numeric>> {
    static void to_json(json& j, const data_view_t<Numeric>& d) {
      for(uint_t ii = 0; ii < uint_t(d.size()); ++ii)
        j.emplace_back(d(ii));
    }
  };

  template <typename Numeric>
  struct adl_serializer<data_span_t<Numeric>> {
    static void to_json(json& j, const data_span_t<Numeric>& d) {
      for(uint_t ii = 0; ii < uint_t(d.size()); ++ii)
        j.emplace_back(d(ii));
    }
  };

  template <typename Numeric>
  struct adl_serializer<data_t<Numeric>> {
    static void from_json(const json& j, data_t<Numeric>& d) {
      d.resize(j.size());
      for(uint_t ii = 0; ii < uint_t(d.size()); ++ii)
        d(ii) = j[ii].get<Numeric>();
    }

    static void to_json(json& j, const data_t<Numeric>& d) {
      for(uint_t ii = 0; ii < uint_t(d.size()); ++ii)
        j.emplace_back(d(ii));
    }
  };

  template <typename Numeric, int N>
  struct adl_serializer<vector_t<Numeric, N>> {
    static void from_json(const json& j, vector_t<Numeric, N>& v) {
      if constexpr(N == n_dynamic)
        v.resize(j.size());

      uint_t len = std::min<uint_t>(v.size(), j.size());
      for(uint_t ii = 0; ii < len; ++ii)
        v(ii) = j[ii].get<Numeric>();
    }

    static void to_json(json& j, const vector_t<Numeric, N>& v) {
      for(uint_t ii = 0; ii < uint_t(v.size()); ++ii)
        j.emplace_back(v(ii));
    }
  };
}    // namespace nlohmann

#endif    //XRD_TYPES_JSON_HPP
