#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestVdShocSolve {

  static void run_property()
  {

  }

  static void run_bfb()
  {
    VdShocSolveData SDS_f90[] = {
      //               shcol, nlev
      VdShocSolveData(1, 4)
    };

    static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(VdShocSolveData);

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize();

      d.cc[0] = 0;
      d.ca[3] = 0;
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    VdShocSolveData SDS_cxx[] = {
      VdShocSolveData(SDS_f90[0]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      vd_shoc_solve(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      vd_shoc_solve_f(d.shcol(), d.nlev(), d.ca, d.cc, d.denom, d.ze, d.var);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    std::cout << std::endl << std::endl << "OUTPUT: " << std::endl;
    for (Int i = 0; i < num_runs; ++i) {
      VdShocSolveData& d_f90 = SDS_f90[i];
      VdShocSolveData& d_cxx = SDS_cxx[i];
      for (Int k = 0; k < d_f90.total1x2(); ++k) {
        //REQUIRE(d_f90.var[k] == d_cxx.var[k]);
        std::cout << d_f90.var[k] << "        " << d_cxx.var[k] << std::endl;
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("vd_shoc_solve_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdShocSolve;

  TestStruct::run_property();
}

TEST_CASE("vd_shoc_solve_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdShocSolve;

  TestStruct::run_bfb();
}

} // namespace
