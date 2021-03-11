#ifndef XRD_ERROR_HPP
#define XRD_ERROR_HPP

#include <system_error>

namespace gsl {
  enum class ec {
    again = -2,
    failure = -1,
    not_in_domain = 1,
    out_of_range,
    invalid_pointer,
    invalid_argument,
    generic_failure,
    factorization_failure,
    sanity_failure,
    no_memory,
    invalid_user_function,
    runaway_iteration,
    max_iterations,
    divide_by_zero,
    invalid_tolerance,
    tolerance_failure,
    underflow,
    overflow,
    loss_of_accuracy,
    roundoff_failure,
    invalid_length,
    matrix_not_square,
    singularity,
    divergence,
    unsupported,
    unimplemented,
    cache_limit_exceeded,
    table_limit_exceeded,
    no_progress,
    no_progress_jacobian,
    tolerance_in_F,
    tolerance_in_X,
    tolerance_in_gradient,
    end_of_file
  };

  namespace details {
    class gsl_category final : public std::error_category {
     public:
      [[nodiscard]] std::string message(int ev) const final;
      [[nodiscard]] const char* name() const noexcept final;
    };
  }    // namespace details

  const details::gsl_category& error_category() noexcept;

  std::error_condition make_error_condition(ec e) noexcept;
}    // namespace gsl

#endif    //XRD_ERROR_HPP
