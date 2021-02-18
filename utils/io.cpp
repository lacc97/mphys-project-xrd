#include "io.hpp"

#include <charconv>
#include <fstream>

#include <fmt/format.h>

#include <nlohmann/json.hpp>

#include "string.hpp"

std::pair<stl::vector<real_t>, stl::vector<real_t>> io::load_csv(std::string_view file) {
  std::pair<stl::vector<real_t>, stl::vector<real_t>> data;

  std::string buf(file);
  std::ifstream fis(buf);
  while(std::getline(fis, buf)) {
    auto line = string::view::trim(buf);
    if(line.starts_with('#'))
      continue;

    auto tokens = string::view::tokenize(line, " \t,;");
    if(tokens.size() < 2)
      return {};

    data.first.push_back(std::strtold(tokens[0].data(), nullptr));
    data.second.push_back(std::strtold(tokens[1].data(), nullptr));
  }

  return data;
}

nlohmann::json io::load_json(std::string_view file) {
  nlohmann::json j;

  std::ifstream fis(std::string{file});
  fis >> j;

  return j;
}

void io::write_csv(std::string_view file, std::tuple<std::span<const real_t>, std::span<const real_t>> data, std::string_view sep) {
  auto x = std::get<0>(data), y = std::get<1>(data);
  if(x.size() != y.size())
    throw std::runtime_error(fmt::format("shape mismatch in data: size of x ({}) != size of y ({})", x.size(), y.size()));

  auto* f = fopen(std::string{file}.c_str(), "w");
  for(uint_t ii = 0; ii < x.size(); ++ii)
    fmt::print(f, "{0}{1}{2}\n", x[ii], sep, y[ii]);
  fclose(f);
}

void io::write_json(std::string_view file, const nlohmann::json& j) {
  std::ofstream fos(std::string{file});
  fos << j;
}
