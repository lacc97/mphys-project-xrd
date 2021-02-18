#include "timer.hpp"

#include <fmt/format.h>
#include <fmt/ostream.h>

hr_timer::hr_timer(time_unit t, std::string_view name) noexcept : m_Unit{t}, m_Name{name}, m_Start{}, m_Elapsed{}, m_Started{false}, m_Paused{false} {}

hr_timer::hr_timer(std::string_view name, time_unit t) noexcept : hr_timer(t, name) {}

void hr_timer::start() noexcept {
  if(!m_Started || m_Paused) {
    m_Elapsed = duration_type::zero();
    m_Started = true;
    m_Paused = false;
    m_Start = std::chrono::high_resolution_clock::now();
  }
}

hr_timer::real_type hr_timer::pause() noexcept {
  if(!m_Paused && m_Started) {
    time_point_type pause_time = std::chrono::high_resolution_clock::now();
    m_Elapsed += std::chrono::duration_cast<duration_type>(pause_time - m_Start);
    m_Paused = true;
  }

  return m_Elapsed.count() * m_Unit.units_per_second;
}

hr_timer::real_type hr_timer::stop() noexcept {
  if(m_Started) {
    if(!m_Paused) {
      time_point_type stop_time = std::chrono::high_resolution_clock::now();
      m_Elapsed += std::chrono::duration_cast<duration_type>(stop_time - m_Start);
    }

    m_Started = false;
  }

  return m_Elapsed.count() * m_Unit.units_per_second;
}

void hr_timer::report(FILE* stream) {
  bool must_continue = false;
  if(m_Started) {
    pause();
    must_continue = true;
  }

  if(m_Name.empty())
    fmt::print(stream, "[{0} {1}]\n", m_Elapsed.count() * m_Unit.units_per_second, m_Unit.unit_display_string);
  else
    fmt::print(stream, "{0} [{1} {2}]\n", m_Name, m_Elapsed.count() * m_Unit.units_per_second, m_Unit.unit_display_string);

  if(must_continue)
    start();
}

void hr_timer::report(std::ostream& os) {
  bool must_continue = false;
  if(m_Started) {
    pause();
    must_continue = true;
  }

  if(m_Name.empty())
    fmt::print(os, "[{0} {1}]\n", m_Elapsed.count() * m_Unit.units_per_second, m_Unit.unit_display_string);
  else
    fmt::print(os, "{0} [{1} {2}]\n", m_Name, m_Elapsed.count() * m_Unit.units_per_second, m_Unit.unit_display_string);

  if(must_continue)
    start();
}

void hr_timer::reset() noexcept {
  m_Elapsed = duration_type::zero();
  m_Started = false;
  m_Paused = false;
}
