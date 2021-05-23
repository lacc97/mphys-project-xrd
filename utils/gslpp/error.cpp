#include "error.hpp"

#include <gsl/gsl_errno.h>

std::string gsl::details::gsl_category::message(int ev) const {
  return gsl_strerror(ev);
}

const char* gsl::details::gsl_category::name() const noexcept {
  return "gsl";
}

const gsl::details::gsl_category& gsl::error_category() noexcept {
  static details::gsl_category g_category;
  return g_category;
}

std::error_code gsl::make_error_code(gsl::ec e) noexcept {
  return std::error_code{static_cast<int>(e), error_category()};
}

std::error_condition gsl::make_error_condition(gsl::ec e) noexcept {
  return std::error_condition{static_cast<int>(e), error_category()};
}
