#ifndef XRD_ARRAY_HPP
#define XRD_ARRAY_HPP

#include <array>
#include <type_traits>

namespace utils {
  namespace details {
    template <typename T, typename... Args, size_t... Is>
    constexpr std::array<T, sizeof...(Is)> make_array(std::index_sequence<Is...>, Args&&... args) noexcept(std::is_nothrow_constructible_v<T, Args...>) {
      return {{(static_cast<void>(Is), T{std::forward<Args>(args)...})...}};
    }
  }

  template <typename T, std::size_t N, typename... Args>
  constexpr std::array<T, N> make_array(Args&&... args) noexcept(std::is_nothrow_constructible_v<T, Args...>) {
    return details::make_array<T>(std::make_index_sequence<N>(), std::forward<Args>(args)...);
  }
}

#endif    //XRD_ARRAY_HPP
