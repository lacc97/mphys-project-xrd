#include "form_factor.hpp"

#include <array>

namespace {
  constexpr std::array<std::array<double, 101>, 25> k_F0 = {
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 26, 0, 0,  0,  0,  0,  0,  0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0,  0,  0,  0,  0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 96, 97, 98, 99, 100},
     {0.005, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.983, 0, 0,      0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 77.972, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 95.956, 96.956, 97.962, 98.962, 99.962},
     {0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.961, 0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 77.937, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 95.895, 96.895, 97.902, 98.903, 99.903},
     {0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.932, 0, 0,      0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 77.892, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 95.809, 96.812, 97.815, 98.817, 99.819},
     {0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.88, 0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 77.809, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 95.664, 96.668, 97.673, 98.678, 99.682},
     {0.025, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.814, 0, 0,     0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,     0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 77.702, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 95.48, 96.487, 97.495, 98.502, 99.509},
     {0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.735, 0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 77.574, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 95.261, 96.271, 97.283, 98.293, 99.302},
     {0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.538, 0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 77.253, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 94.726, 95.743, 96.769, 97.784, 98.799},
     {0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.302, 0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 76.865, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 94.109, 95.133, 96.174, 97.196, 98.217},
     {0.07, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24.707, 0, 0,      0,      0,      0,     0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,     0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 75.854, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 92.573, 93.611, 94.715, 95.75, 96.784},
     {0.09, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24.023, 0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 74.638, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 90.886, 91.934, 93.116, 94.162, 95.207},
     {0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 23.666, 0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 73.976, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 89.998, 91.051, 92.271, 93.322, 94.372},
     {0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 22.736, 0, 0,      0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 72.162, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 87.761, 88.819, 90.134, 91.193, 92.252},
     {0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21.803, 0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 70.212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 85.442, 86.499, 87.881, 88.943, 90.007},
     {0.175, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20.911, 0, 0,      0,     0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,     0,      0,      0,     0,
      0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 68.22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 83.117, 84.17, 85.563, 86.627, 87.694},
     {0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20.033, 0, 0,      0,      0,      0,     0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,     0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 66.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 80.766, 81.807, 83.182, 84.24, 85.304},
     {0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18.34, 0, 0,      0,      0,      0,      0,    0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0,      0,      0,      0,      0,    0,
      0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 62.139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 76.041, 77.038, 78.314, 79.343, 80.38},
     {0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16.729, 0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 58.291, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 71.555, 72.488, 73.624, 74.598, 75.585},
     {0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13.809, 0, 0,     0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,     0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 51.466, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 63.69, 64.458, 65.315, 66.134, 66.973},
     {0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11.468, 0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 45.893, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 57.217, 57.851, 58.491, 59.165, 59.859},
     {0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.7165, 0, 0,      0,      0,      0,      0,      0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,      0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 41.237, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 51.701, 52.248, 52.753, 53.323, 53.9004},
     {0.7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.4697, 0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 37.192, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 46.943, 47.433, 47.867, 48.363, 48.865},
     {0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.6042, 0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 33.589, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 42.878, 43.331, 43.737, 44.182, 44.627},
     {0.9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.9889, 0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0,      0,      0,      0,      0,     0,
      0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 30.347, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 39.421, 39.851, 40.253, 40.667, 41.076},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.515, 0, 0,     0,      0,      0,      0,     0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0,     0,      0,      0,      0,     0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27.439, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 36.44, 36.862, 37.272, 37.673, 38.065}}};
  constexpr size_t k_F0_element_count = std::tuple_size_v<std::tuple_element_t<0, decltype(k_F0)>> - 1;

  const auto g_f0_splines = []() -> std::array<gsl::spline, k_F0_element_count> {
    auto fn_column = [](size_t c) constexpr noexcept->std::array<double, k_F0.size()> {
      std::array<double, k_F0.size()> arr;
      for(size_t ii = 0; ii < arr.size(); ++ii)
        arr[ii] = k_F0[ii][c];
      return arr;
    };

    auto x = fn_column(0);

    return [&x, &fn_column]<size_t... Is>(std::index_sequence<Is...>)->std::array<gsl::spline, sizeof...(Is)> {
      return {{gsl::spline(gsl::interp_type::steffen, x, fn_column(Is + 1))...}};
    }
    (std::make_index_sequence<k_F0_element_count>());
  }();
  thread_local auto g_f0_accels = std::array<gsl::interp_accel, k_F0_element_count>{};
}    // namespace

//#include <fmt/format.h>

real_t xrd::tables::f0(uint_t Z, real_t x) noexcept {
//  return Z;
//  fmt::print("f_0({}, {}) -> {}\n", Z, x, g_f0_splines[Z - 1](x, g_f0_accels[Z - 1]));
  return g_f0_splines[Z - 1](x, g_f0_accels[Z - 1]);
}
