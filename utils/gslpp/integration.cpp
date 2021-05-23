#include "integration.hpp"

#include <gsl/gsl_integration.h>

#include "error.hpp"

namespace {
  int from_rule(gsl::integration::rule r) noexcept {
    switch(r) {
      case gsl::integration::rule::e_Gauss_15:
        return GSL_INTEG_GAUSS15;
      case gsl::integration::rule::e_Gauss_21:
        return GSL_INTEG_GAUSS21;
      default:
      case gsl::integration::rule::e_Gauss_31:
        return GSL_INTEG_GAUSS31;
      case gsl::integration::rule::e_Gauss_41:
        return GSL_INTEG_GAUSS41;
      case gsl::integration::rule::e_Gauss_51:
        return GSL_INTEG_GAUSS51;
      case gsl::integration::rule::e_Gauss_61:
        return GSL_INTEG_GAUSS61;
    }
  }
}    // namespace

namespace gsl::integration {
  namespace {
    class workspace {
      struct destructor {
        inline void operator()(gsl_integration_workspace* w) const noexcept {
          if(w)
            gsl_integration_workspace_free(w);
        }
      };

     public:
      explicit workspace(size_t n) : m_N{n}, mp_ptr{gsl_integration_workspace_alloc(m_N)} {
        if(!mp_ptr)
          throw std::system_error{make_error_code(ec::failure)};
      }

      void resize(size_t n) {
        if(m_N < n) {
          mp_ptr.reset(gsl_integration_workspace_alloc(n));
          if(!mp_ptr)
            throw std::system_error{make_error_code(ec::failure)};
          m_N = n;
        }
      }

      operator gsl_integration_workspace*() noexcept {
        return mp_ptr.get();
      }

     private:
      size_t m_N;
      std::unique_ptr<gsl_integration_workspace, destructor> mp_ptr;
    };
  }    // namespace
}    // namespace gsl::integration


gsl::integration::result gsl::integration::detail::qag_adaptive(gsl::integration::detail::function f, double a, double b, gsl::integration::rule r,
                                                                             double epsabs, double epsrel, size_t limit) {
  gsl::integration::workspace workspace{limit == 0 ? 64 : limit};

  gsl_function fun{f.function, const_cast<void*>(f.params)};

  result res{};
  if(auto status = gsl_integration_qag(&fun, a, b, epsabs, epsrel, limit, from_rule(r), workspace, &res.res, &res.abserr); status != 0)
    throw std::system_error{status, error_category()};
  return res;
}
