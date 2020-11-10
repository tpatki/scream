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

  static void run_decomp()
  {

  }

  static void run_solver()
  {
//    const bool random = true;
//    static constexpr Real var[3] = {10, 5, 15};
//    static constexpr Real kv_term[3] = {2, 6, 4};
//    static constexpr Real tmpi[3] = {1, 5, 1};
//    static constexpr Real rdp_zt[3] = {15, 3, 2};
//    static constexpr Real flux[2] = {1, 4};

//    static constexpr int shcol = 10;
//    static constexpr int nlev = 71;

//    VdShocSolveData  f90_solve(shcol, nlev, 1);
//    VdShocDecompData f90_decomp(shcol,nlev,nlev+1, 1);
//    VdShocSolveData cxx_solve(shcol,nlev,1);

//    if (random) {
//      f90_solve.randomize({ {f90_solve.var, {1,2}} });
//      f90_decomp.randomize({{f90_decomp.kv_term, {1, 2}},
//                            {f90_decomp.tmpi, {1,2}},
//                            {f90_decomp.rdp_zt,{1,2}},
//                            {f90_decomp.flux, {1,2}} });
//    }
//    else {
//      // Fill in test data
//      for(Int s = 0; s < shcol; ++s) {
//        // First on the nlev grid
//        for(Int n = 0; n < nlev; ++n) {
//          const auto index = n + s * nlev;
//          f90_solve.var[index] = var[n];
//          cxx_solve.var[index] = var[n];

//          f90_decomp.kv_term[index] = kv_term[n];
//          f90_decomp.tmpi[index] = tmpi[n];
//          f90_decomp.rdp_zt[index] = rdp_zt[n];
//        }
//        f90_decomp.flux[s] = flux[s];
//      }
//    }

//    vd_shoc_decomp(f90_decomp);

//    // Fill in test data
//    for(Int s = 0; s < shcol; ++s) {
//      // First on the nlev grid
//      for(Int n = 0; n < nlev; ++n) {
//        const auto index = n + s * nlev;
//        f90_solve.ca[index] = f90_decomp.ca[index];
//        f90_solve.cc[index] = f90_decomp.cc[index];
//        f90_solve.denom[index] = f90_decomp.denom[index];
//        f90_solve.ze[index] = f90_decomp.ze[index];
//        f90_solve.rdp_zt[index] = f90_decomp.rdp_zt[index];

//        cxx_solve.ca[index] = f90_decomp.ca[index];
//        cxx_solve.cc[index] = f90_decomp.cc[index];
//        cxx_solve.denom[index] = f90_decomp.denom[index];
//        cxx_solve.ze[index] = f90_decomp.ze[index];
//        cxx_solve.rdp_zt[index] = f90_decomp.rdp_zt[index];
//        cxx_solve.var[index] = f90_solve.var[index];
//      }

//      f90_solve.flux[s] = f90_decomp.flux[s];
//      cxx_solve.flux[s] = f90_decomp.flux[s];
//    }

//    vd_shoc_solve(f90_solve);

//    cxx_solve.transpose<ekat::TransposeDirection::c2f>();
//    vd_shoc_solve_f(cxx_solve.shcol(), cxx_solve.nlev(), cxx_solve.ca, cxx_solve.cc, cxx_solve.denom, cxx_solve.ze,
//                    cxx_solve.rdp_zt, cxx_solve.dtime, cxx_solve.flux,
//                    cxx_solve.var);
//    cxx_solve.transpose<ekat::TransposeDirection::f2c>();

//    // Verify BFB results, all data should be in C layout
//    std::cout << std::endl << std::endl << "OUTPUT: " << std::endl;
//    for (Int k = 0; k < f90_solve.total1x2(); ++k) {
//      std::cout << f90_solve.var[k] << "        " << cxx_solve.var[k] << std::endl;
//      const Real rel_diff
//          = std::abs(f90_solve.var[k] - cxx_solve.var[k])/(std::abs(f90_solve.var[k]) == 0 ? 1 : std::abs(f90_solve.var[k]));
//      REQUIRE(rel_diff < 1e-5);
//    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("vd_shoc_solve_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdShocSolve;

  TestStruct::run_decomp();
}

TEST_CASE("vd_shoc_solve_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdShocSolve;

  TestStruct::run_solver();
}

} // namespace
