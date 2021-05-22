#ifndef XRD_IO_HPP
#define XRD_IO_HPP

#include <span>

#include <nlohmann/json_fwd.hpp>

#include "data/dataset_2d.hpp"
#include "templates/array.hpp"
#include "types.hpp"

namespace io {
  namespace detail {
    void write_csv(std::string_view file, std::span<const std::span<const real_t>> datasets, std::string_view sep = " ");
  }

  std::pair<stl::vector<real_t>, stl::vector<real_t>> load_csv(std::string_view file);

  void write_csv(std::string_view file, std::tuple<std::span<const real_t>, std::span<const real_t>> data, std::string_view sep = " ");
  template <typename T1, typename T2>
  inline void write_csv(std::string_view file, std::tuple<T1, T2> data, std::string_view sep = " ") {
    return write_csv(file, std::make_tuple<std::span<const real_t>, std::span<const real_t>>(std::get<0>(data), std::get<1>(data)), sep);
  }
  inline void write_csv(std::string_view file, ds::dataset_2d_view data, std::string_view sep = " ") {
    return write_csv(file, std::make_tuple<std::span<const real_t>, std::span<const real_t>>(data.x(), data.y()), sep);
  }

  template <typename... Args>
  void write_csv(std::string_view file, Args&&... args) {
    detail::write_csv(file, utils::make_array_of<std::span<const real_t>>(std::forward<Args>(args)...));
  }


  nlohmann::json load_json(std::string_view file);
  void write_json(std::string_view file, const nlohmann::json& j);
}    // namespace io

#endif    //XRD_IO_HPP
