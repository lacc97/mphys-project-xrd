#ifndef XRD_ALLOCATOR_HPP
#define XRD_ALLOCATOR_HPP

#include <memory>

namespace alloc {
  // Allocator adaptor that interposes construct() calls to
  // convert value initialization into default initialization.
  // From: https://stackoverflow.com/questions/21028299/is-this-behavior-of-vectorresizesize-type-n-under-c11-and-boost-container/21028912#21028912
  template <typename T, typename A = std::allocator<T>>
  class default_init_allocator : public A {
    typedef std::allocator_traits<A> a_t;

   public:
    template <typename U>
    struct rebind {
      using other = default_init_allocator<U, typename a_t::template rebind_alloc<U>>;
    };

    using A::A;

    template <typename U>
    inline void construct(U* ptr) noexcept(std::is_nothrow_default_constructible_v<U>) {
      ::new(static_cast<void*>(ptr)) U;
    }

    template <typename U, typename... Args>
    inline void construct(U* ptr, Args&&... args) {
      a_t::construct(static_cast<A&>(*this), ptr, std::forward<Args>(args)...);
    }
  };
}    // namespace alloc

#endif    //XRD_ALLOCATOR_HPP
