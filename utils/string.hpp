#ifndef XRD_STRING_HPP
#define XRD_STRING_HPP

#include "types.hpp"

namespace string {
  namespace view {
    stl::vector<std::string_view> tokenize(std::string_view s, std::string_view delims);

    std::string_view ltrim(std::string_view s, std::string_view chars = "\t\n\v\f\r ") noexcept;
    std::string_view rtrim(std::string_view s, std::string_view chars = "\t\n\v\f\r ") noexcept;
    inline std::string_view trim(std::string_view s, std::string_view chars = "\t\n\v\f\r ") noexcept {
      return rtrim(ltrim(s, chars), chars);
    }
  }    // namespace view

  inline stl::vector<std::string> tokenize(std::string_view s, std::string_view delims) {
    auto token_views = view::tokenize(s, delims);
    return stl::vector<std::string>(token_views.begin(), token_views.end());
  }

  inline std::string ltrim(std::string_view s, std::string_view chars = "\t\b\v\f\r ") {
    return std::string(view::ltrim(s, chars));
  }
  inline std::string rtrim(std::string_view s, std::string_view chars = "\t\b\v\f\r ") {
    return std::string(view::rtrim(s, chars));
  }
  inline std::string trim(std::string_view s, std::string_view chars = "\t\b\v\f\r ") {
    return std::string(view::trim(s, chars));
  }
}    // namespace string

#endif    //XRD_STRING_HPP
