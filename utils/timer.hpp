#pragma once

#include <chrono>
#include <iostream>
#include <string>

struct time_unit {
  long double units_per_second;
  std::string_view unit_display_string;
};

constexpr time_unit k_Hours{1.0 / 3600.0, "h"};
constexpr time_unit k_Minutes{1.0 / 60.0, "m"};
constexpr time_unit k_Seconds{1, "s"};
constexpr time_unit k_Milliseconds{1000, "ms"};
constexpr time_unit k_Microseconds{1000000, "μs"};
constexpr time_unit k_Nanoseconds{1000000000L, "ns"};    // likely inaccurate?

class hr_timer {
  using real_type = long double;
  using duration_type = std::chrono::duration<real_type>;
  using time_point_type = std::chrono::high_resolution_clock::time_point;

 public:
  explicit hr_timer(time_unit t = k_Milliseconds, std::string_view name = "") noexcept;
  explicit hr_timer(std::string_view name, time_unit t = k_Milliseconds) noexcept;

  void start() noexcept;
  real_type pause() noexcept;
  real_type stop() noexcept;

  void reset() noexcept;

  void report(FILE* stream = stdout);
  void report(std::ostream& os);

 private:
  time_unit m_Unit;
  std::string m_Name;
  time_point_type m_Start;
  duration_type m_Elapsed;
  bool m_Started;
  bool m_Paused;
};

template <typename F>
void benchmark(std::string_view name, F&& f) {
  hr_timer timer{name};
  timer.start();
  f();
  timer.stop();
  timer.report();
}

template <typename F>
void benchmark(std::string_view name, time_unit t, F&& f) {
  hr_timer timer{name, t};
  timer.start();
  f();
  timer.stop();
  timer.report();
}
