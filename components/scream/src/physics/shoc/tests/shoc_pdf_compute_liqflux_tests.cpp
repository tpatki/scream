#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_arch.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocPdfCompLiqFlux {

  static void run_property()
  {
  
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_liquid_water_flux

    // TEST ONE
    // Dry test.  Verify that wqls is zero if there is no
    //  liquid water loading.
    
    // Input data (note that many of these inputs will be used
    //  for test two as well)
    
    // Define the gaussian fraction [-]
    static constexpr Real a = 0.5;
    // Define vertical velocity of first gaussian
    static constexpr Real w1_1 = 1;
    // Define liquid water loading gaussian one
    static constexpr Real ql1 = 0;
    // Define vertical velocity of second gaussian
    static constexpr Real w1_2 = -1;
    // Define liquid water loading gaussian two
    static constexpr Real ql2 = 0;
  
    // Initialize data structure for bridging to F90
    SHOCPDFcompliqfluxData SDS;
    
    // Load input data
    SDS.a = a;
    SDS.w1_1 = w1_1;
    SDS.ql1 = ql1;
    SDS.w1_2 = w1_2;
    SDS.ql2 = ql1;
    
    // Verify there is no liquid water in inputs
    REQUIRE(SDS.ql1 == 0);
    REQUIRE(SDS.ql2 == 0);
    
    // Call the fortran implementation
    shoc_assumed_pdf_compute_liquid_water_flux(SDS);
    
    // Verify that liquid water flux is zero
    REQUIRE(SDS.wqls == 0);
    
    // TEST TWO
    // Convective test.  Given inputs representative of convective
    //  conditions, verify that liquid water flux is greater than zero.
    
     // Define the gaussian fraction [-]
    static constexpr Real a_conv = 0.1;
    // Define vertical velocity of first gaussian
    static constexpr Real w1_1_conv = 3;
    // Define liquid water loading gaussian one
    static constexpr Real ql1_conv = 0.004;
    // Define vertical velocity of second gaussian
    static constexpr Real w1_2_conv = -0.1;
    // Define liquid water loading gaussian two
    static constexpr Real ql2_conv = 0; 
    
    // Load input data
    SDS.a = a_conv;
    SDS.w1_1 = w1_1_conv;
    SDS.ql1 = ql1_conv;
    SDS.w1_2 = w1_2_conv;
    SDS.ql2 = ql1_conv;  
    
    // Verify input data is as we expect for convective conditions
    REQUIRE(SDS.a < 0.5);
    REQUIRE(SDS.w1_1 > 0);
    REQUIRE(SDS.w1_2 < 0);
    
    // Call the fortran implementation
    shoc_assumed_pdf_compute_liquid_water_flux(SDS);
    
    // Verify the result.  Liquid water flux should be positive
    REQUIRE(SDS.wqls > 0);
    
    // TEST THREE
    // negatively skewed test.  Given conditions indicative of negative
    //  skewness, with condensate in both gaussians, verify that liquid
    //  water flux is negative
    
     // Define the gaussian fraction [-]
    static constexpr Real a_neg = 0.7;
    // Define vertical velocity of first gaussian
    static constexpr Real w1_1_neg = 0.1;
    // Define liquid water loading gaussian one
    static constexpr Real ql1_neg = 0.004;
    // Define vertical velocity of second gaussian
    static constexpr Real w1_2_neg = -2.0;
    // Define liquid water loading gaussian two
    static constexpr Real ql2_neg = 0.004; 
    
    // Load input data
    SDS.a = a_neg;
    SDS.w1_1 = w1_1_neg;
    SDS.ql1 = ql1_neg;
    SDS.w1_2 = w1_2_neg;
    SDS.ql2 = ql2_neg;  
    
    // Verify input data is as we expect for convective conditions
    REQUIRE(SDS.a > 0.5);
    REQUIRE(SDS.w1_1 > 0);
    REQUIRE(SDS.w1_2 < 0);
    REQUIRE(SDS.ql1 == SDS.ql2);
    
    // Call the fortran implementation
    shoc_assumed_pdf_compute_liquid_water_flux(SDS);
    
    // Verify the result.  Liquid water flux should be positive
    REQUIRE(SDS.wqls < 0);
        
  }
  
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_pdf_compute_liqflux_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfCompLiqFlux;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_compute_liqflux_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfCompLiqFlux;

}

} // namespace
