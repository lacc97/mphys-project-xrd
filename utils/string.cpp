#include "string.hpp"

stl::vector<std::string_view> string::view::tokenize(std::string_view s, std::string_view delims) {
  stl::vector<std::string_view> tokens;

  while(!s.empty()) {
    auto delimPos = s.find_first_of(delims);
    if(delimPos == std::string_view::npos) {
      tokens.push_back(s);
      break;
    } else {
      auto token = s.substr(0, delimPos);
      if(!token.empty())
        tokens.push_back(token);
      s = s.substr(delimPos + 1);
    }
  }

  return tokens;
}
std::string_view string::view::ltrim(std::string_view s, std::string_view chars) noexcept {
  auto pos = s.find_first_not_of(chars);
  if(pos == std::string_view::npos)
    return {};
  return s.substr(pos);
}
std::string_view string::view::rtrim(std::string_view s, std::string_view chars) noexcept {
  auto pos = s.find_last_not_of(chars);
  if(pos == std::string_view::npos)
    return {};
  return s.substr(0, pos + 1);
}
