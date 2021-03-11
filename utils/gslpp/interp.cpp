#include "interp.hpp"

#include <gsl/gsl_interp.h>

#include "error.hpp"

void gsl::interp_accel::interp_accel_deleter::operator()(void* ptr) const noexcept {
  if(ptr)
    gsl_interp_accel_free(static_cast<gsl_interp_accel*>(ptr));
}

gsl::interp_accel::interp_accel() : mp_accel{gsl_interp_accel_alloc()} {}

size_t gsl::interp_accel::find(std::span<double> array, double x) noexcept {
  return gsl_interp_accel_find(static_cast<gsl_interp_accel*>(mp_accel.get()), array.data(), array.size(), x);
}

void gsl::interp_accel::reset() {
  if(auto status = gsl_interp_accel_reset(static_cast<gsl_interp_accel*>(mp_accel.get())); status)
    throw std::system_error{status, error_category()};
}
