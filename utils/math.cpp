#include "math.hpp"

void math::data::details::linspace(real_t start, real_t stop, std::span<real_t> out) {
  if(out.size() <= 1)
    throw std::invalid_argument(fmt::format("invalid count ({}): must be at least 2", out.size()));
  if(stop <= start)
    throw std::invalid_argument(fmt::format("invalid endpoints: start ({}) must be smaller than stop ({})", start, stop));

  const real_t h = (stop - start) / static_cast<real_t>(out.size() - 1);

#pragma omp parallel for
  for(size_t ii = 0; ii < out.size() - 1; ++ii)
    out[ii] = start + ii * h;
  out.back() = stop;
}

thread_local std::mt19937_64 math::rand::tl_Generator{std::random_device{}()};

real_t math::rand::unit() {
  thread_local std::uniform_real_distribution<real_t> tl_Uniform(0, 1);
  return tl_Uniform(tl_Generator);
}

rmatrix_t<3, 3> math::linalg::generate_orthonormal_basis_from_vector(const rvec3_t& z_axis) noexcept {
  rvec3_t x_axis = z_axis.cross(rvec3_t{1, 0, 0});
  if(x_axis.squaredNorm() == 0)
    x_axis = z_axis.cross(rvec3_t{0, 1, 0});
  assert(x_axis.squaredNorm() != 0);

  rvec3_t y_axis = z_axis.cross(x_axis);
  assert(y_axis.squaredNorm() != 0);

  rmatrix_t<3, 3> basis_mat;
  basis_mat.col(0) = x_axis.normalized();
  basis_mat.col(1) = y_axis.normalized();
  basis_mat.col(2) = z_axis.normalized();

  return basis_mat;
}
