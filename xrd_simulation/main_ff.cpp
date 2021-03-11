#include <io.hpp>
#include <math.hpp>
#include <types.hpp>

#include "tables/form_factor.hpp"

int main(int argc, char** argv) {
  auto fn_gen_ff = [](uint_t Z, const rdata_t& x) {
    rdata_t ff{x.size()};
    for(size_t ii = 0; ii < ff.size(); ++ii)
      ff(ii) = xrd::tables::f0(Z, x(ii));
    return ff;
  };

  const rdata_t x = math::data::linspace(0, 1, 100000);
  const rdata_t fe = fn_gen_ff(26, x);
  const rdata_t pt = fn_gen_ff(78, x);

  io::write_csv("fept/fe_form_factor.csv", std::tie(x, fe));
  io::write_csv("fept/pt_form_factor.csv", std::tie(x, pt));
}