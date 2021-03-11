#ifndef XRD_COMPLEX_HPP
#define XRD_COMPLEX_HPP

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace math {
  class complex {
   public:
    constexpr complex(double r = 0.0) noexcept : complex{r, 0} {}
    constexpr complex(double r, double i) noexcept : m_value{r, i} {}
    constexpr complex(gsl_complex z) noexcept : m_value{z} {}

    [[nodiscard]] constexpr double re() const noexcept {
      return m_value.dat[0];
    }
    [[nodiscard]] constexpr double im() const noexcept {
      return m_value.dat[1];
    }

    [[nodiscard]] constexpr complex conj() const noexcept {
      return {re(), -im()};
    }
    [[nodiscard]] constexpr double norm() const noexcept {
      return std::sqrt(squared_norm());
    }
    [[nodiscard]] constexpr double squared_norm() const noexcept {
      return re() * re() + im() * im();
    }

    constexpr complex operator-() const noexcept {
      return {-re(), -im()};
    }
    constexpr complex operator+(complex z) const noexcept {
      return {re() + z.re(), im() + z.im()};
    }
    constexpr complex operator+(double r) const noexcept {
      return {re() + r, im()};
    }
    constexpr complex& operator+=(complex z) noexcept {
      m_value = {re() + z.re(), im() + z.im()};
      return *this;
    }
    constexpr complex& operator+=(double z) noexcept {
      m_value = {re() + z, im()};
      return *this;
    }
    constexpr complex operator-(complex z) const noexcept {
      return {re() - z.re(), im() - z.im()};
    }
    constexpr complex& operator-=(complex z) noexcept {
      m_value = {re() - z.re(), im() - z.im()};
      return *this;
    }
    constexpr complex& operator-=(double z) noexcept {
      m_value = {re() - z, im()};
      return *this;
    }
    constexpr complex operator*(complex z) const noexcept {
      return {re() * z.re() - im() * z.im(), re() * z.im() + im() * z.re()};
    }
    constexpr complex operator*(double r) const noexcept {
      return {re() * r, im() * r};
    }
    constexpr complex& operator*=(complex z) noexcept {
      m_value = {re() * z.re() - im() * z.im(), re() * z.im() + im() * z.re()};
      return *this;
    }
    constexpr complex& operator*=(double r) noexcept {
      m_value = {re() * r, im() * r};
      return *this;
    }
    constexpr complex operator/(complex z) const noexcept {
      return complex{re() * z.re() + im() * z.im(), re() * z.im() - im() * z.re()} / squared_norm();
    }
    constexpr complex operator/(double r) const noexcept {
      return {re() / r, im() / r};
    }
    constexpr complex& operator/=(complex z) noexcept {
      m_value = {re() * z.re() + im() * z.im(), re() * z.im() - im() * z.re()};
      return operator/=(squared_norm());
    }
    constexpr complex& operator/=(double r) noexcept {
      m_value = {re() / r, im() / r};
      return *this;
    }

    constexpr explicit operator gsl_complex() const noexcept {
      return m_value;
    }

   private:
    gsl_complex m_value;
  };

  inline constexpr complex operator+(double a, complex b) noexcept {
    return b + a;
  }
  inline constexpr complex operator-(double a, complex b) noexcept {
    return -(b - a);
  }
  inline constexpr complex operator*(double a, complex b) noexcept {
    return b * a;
  }
  inline constexpr complex operator/(double a, complex b) noexcept {
    return complex{a} / b;
  }
}    // namespace math

#endif    //XRD_COMPLEX_HPP
