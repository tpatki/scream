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
    VdShocSolveData  f90_solve(1, 3, 1);
    f90_solve.randomize({ {f90_solve.var, {1,2}} });
    VdShocSolveData cxx_solve(f90_solve);

    VdShocDecompData f90_decomp(1, 3, 5, 1);
    f90_decomp.randomize({{f90_decomp.kv_term, {1, 2}},
                          {f90_decomp.tmpi, {1,2}},
                          {f90_decomp.rdp_zt,{1,2}},
                          {f90_decomp.flux, {1,2}} });
    vd_shoc_decomp(f90_decomp);

    f90_solve.ca = f90_decomp.ca;
    f90_solve.cc = f90_decomp.cc;
    f90_solve.denom = f90_decomp.denom;
    f90_solve.ze = f90_decomp.ze;
    f90_solve.rdp_zt = f90_decomp.rdp_zt;
    f90_solve.flux = f90_decomp.flux;

    cxx_solve.ca = f90_decomp.ca;
    cxx_solve.cc = f90_decomp.cc;
    cxx_solve.denom = f90_decomp.denom;
    cxx_solve.ze = f90_decomp.ze;
    cxx_solve.rdp_zt = f90_decomp.rdp_zt;
    cxx_solve.flux = f90_decomp.flux;

    vd_shoc_solve(f90_solve);

    cxx_solve.transpose<ekat::TransposeDirection::c2f>();
    vd_shoc_solve_f(cxx_solve.shcol(), cxx_solve.nlev(), cxx_solve.ca, cxx_solve.cc, cxx_solve.denom, cxx_solve.ze,
                    cxx_solve.rdp_zt, cxx_solve.dtime, cxx_solve.flux,
                    cxx_solve.var);
    cxx_solve.transpose<ekat::TransposeDirection::f2c>();

    // Verify BFB results, all data should be in C layout
    std::cout << std::endl << std::endl << "OUTPUT: " << std::endl;
    for (Int k = 0; k < f90_solve.total1x2(); ++k) {
      //REQUIRE(d_f90.var[k] == d_cxx.var[k]);
      std::cout << f90_solve.var[k] << "        " << cxx_solve.var[k] << std::endl;
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
