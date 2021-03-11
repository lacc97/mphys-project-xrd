#include "spline.hpp"

#include <stdexcept>

#include <fmt/format.h>

#include <gsl/gsl_spline.h>

#include "error.hpp"
#include "interp_p.hpp"

void gsl::spline::spline_deleter::operator()(void* ptr) const noexcept {
  if(ptr)
    gsl_spline_free(static_cast<gsl_spline*>(ptr));
}

double gsl::spline::d(double x, gsl::interp_accel& acc) const noexcept {
  return gsl_spline_eval_deriv(static_cast<gsl_spline*>(mp_spline.get()), x, static_cast<gsl_interp_accel*>(acc.mp_accel.get()));
}

double gsl::spline::d2(double x, gsl::interp_accel& acc) const noexcept {
  return gsl_spline_eval_deriv2(static_cast<gsl_spline*>(mp_spline.get()), x, static_cast<gsl_interp_accel*>(acc.mp_accel.get()));
}

void gsl::spline::init(std::span<const double> x, std::span<const double> y) {
  if(x.size() != y.size())
    throw std::invalid_argument(fmt::format("shape mismatch: x and y have different sizes ({} != {})", x.size(), y.size()));

  mp_spline.reset(gsl_spline_alloc(details::get_type(m_type), x.size()));
  if(auto status = gsl_spline_init(static_cast<gsl_spline*>(mp_spline.get()), x.data(), y.data(), x.size()); status)
    throw std::system_error{status, error_category()};
}

double gsl::spline::operator()(double x, gsl::interp_accel& acc) const noexcept {
  return gsl_spline_eval(static_cast<gsl_spline*>(mp_spline.get()), x, static_cast<gsl_interp_accel*>(acc.mp_accel.get()));
}
