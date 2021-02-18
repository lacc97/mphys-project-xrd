#ifndef XRD_SIMULATED_ANNEALING_HPP
#define XRD_SIMULATED_ANNEALING_HPP

#include "math.hpp"
#include "types.hpp"

#include <fmt/format.h>

namespace opt {
  template <typename AnnealingTraits, bool Trace = false>
  class simulated_annealer {
    using traits_type = AnnealingTraits;
    using solution_type = typename traits_type::solution_type;

   public:
    solution_type run(const traits_type& traits, const sint_t num_iterations, const sint_t steps) const {
      if(num_iterations <= 0)
        throw std::runtime_error(fmt::format("invalid num_iterations ({})", num_iterations));
      if(steps <= 0)
        throw std::runtime_error(fmt::format("invalid steps ({})", steps));

      solution_type sMin = traits.initial_solution();
      long double eMin = traits.energy(sMin);

      if constexpr(Trace)
        fmt::print("Started simulated annealing with {0} iterations and {1} steps per iteration:\n", num_iterations, steps);

      for(sint_t k = 0; k < num_iterations; ++k) {
        solution_type s = traits.initial_solution();
        real_t e = traits.energy(s);

        for(sint_t l = 0; l < steps; ++l) {
          // Calculate the temperature.
          real_t t = T(l);

          // Pick a random neighbouring solution.
          solution_type sNew = traits.random_neighbour(s);
          real_t ep = traits.energy(sNew);

          // Attempt to jump.
          if(P(e, ep, t) >= math::rand::unit()) {
            s = sNew;
            e = ep;

            // Check if the new solution is better.
            if(e < eMin) {
              sMin = s;
              eMin = e;

              if constexpr(Trace) {
                if constexpr(fmt::has_formatter<solution_type, fmt::format_context>::value)
                  fmt::print("  On step {0} of iteration {1}, a better solution {2} with energy {3} was found.\n", l + 1, k + 1, sMin, eMin);
                else
                  fmt::print("  On step {0} of iteration {1}, a better solution with energy {2} was found.\n", l + 1, k + 1, eMin);
              }
            }
          }
        }
      }

      return sMin;
    }

   private:
    /*
    * This static method returns the probability of jumping from a solution  with
    * energy e to a solution with energy ep at temperature t. This probability function
    * is the same as the default used by MATLAB
    * (https://www.mathworks.com/help/gads/how-simulated-annealing-works.html).
    */
    inline static real_t P(real_t e, real_t ep, real_t t) {
      // Always jump if we have found a better solution.
      if(ep < e)
        return 1;

      // Otherwise, explore the rest of the search space with a
      // probability that decreases as the temperature t decreases.
      // This probability lies between 0 and 1/2.
      return 1 / (1 + std::exp((ep - e) / t));
    }

    /*
    * This static method returns the temperature as a function of the iteration
    * number. This temperature function is based on the one use by MATLAB
    * (https://www.mathworks.com/help/gads/how-simulated-annealing-works.html).
    */
    inline static real_t T(sint_t l) {
      return 273.15 * std::pow(0.999, l);
    }
  };
}    // namespace opt

#endif    //XRD_SIMULATED_ANNEALING_HPP
