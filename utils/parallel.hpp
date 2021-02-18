#ifndef COMPUTINGPROJECT_PARALLEL_HPP
#define COMPUTINGPROJECT_PARALLEL_HPP

#include <stdexcept>

namespace Parallel {
  template <typename IntT>
  class IntRange {
    using value_type = IntT;

    class Iterator {};

   public:
    inline IntRange(value_type first, value_type last, value_type step = 1) : m_Begin{first}, m_End{last}, m_Step{step} {
      if(m_Step == 0)
        throw std::invalid_argument("step cannot be 0");
    }

   private:
    value_type m_Begin;
    value_type m_End;
    value_type m_Step;
  };
}    // namespace Parallel

#endif    //COMPUTINGPROJECT_PARALLEL_HPP
