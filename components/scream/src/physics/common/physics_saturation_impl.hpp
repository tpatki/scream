#ifndef PHYSICS_SATURATION_IMPL_HPP
#define PHYSICS_SATURATION_IMPL_HPP

#include "physics_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace physics {

/*
 * Implementation of saturation functions. Clients should NOT #include
 * this file, #include physics_functions.hpp instead.
 */


template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::svp_murphy_koop(const Spack& t, const bool ice)
{
  // Murphy & Koop (2005)
  Spack result;
  const auto tmelt = C::Tmelt;
  Smask ice_mask = (t < tmelt) && ice;
  Smask liq_mask = !ice_mask;

  if (ice_mask.any()) {
    // (good down to 110 K)
    //creating array for storing coefficients of ice sat equation
    const Scalar ic[]= {9.550426, 5723.265, 3.53068, 0.00728332};
    Spack ice_result = exp(ic[0] - (ic[1] / t) + (ic[2] * log(t)) - (ic[3] * t));

    result.set(ice_mask, ice_result);
  }

  if (liq_mask.any()) {
    // (good for 123 < T < 332 K)

    //creating array for storing coefficients of liq sat equation
    const Scalar lq[] = {54.842763, 6763.22, 4.210, 0.000367, 0.0415, 218.8, 53.878,
			 1331.22, 9.44523, 0.014025 };

    Spack liq_result = exp(lq[0] - (lq[1] / t) - (lq[2] * log(t)) + (lq[3] * t) +
			   (tanh(lq[4] * (t - lq[5])) * (lq[6] - (lq[7] / t) -
							 (lq[8] * log(t)) + lq[9] * t)));

    result.set(liq_mask, liq_result);
  }

  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::polysvp1(const Spack& t, const bool ice)
{
  // REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

  // ice
  static Scalar ai[] = {
    6.11147274,     0.503160820,     0.188439774e-1,
    0.420895665e-3, 0.615021634e-5,  0.602588177e-7,
    0.385852041e-9, 0.146898966e-11, 0.252751365e-14};

  // liquid, V1.7
  static Scalar a[] = {
    6.11239921,      0.443987641,     0.142986287e-1,
    0.264847430e-3,  0.302950461e-5,  0.206739458e-7,
    0.640689451e-10,-0.952447341e-13,-0.976195544e-15};

  Spack dt = pack::max(t - sp(273.15), sp(-80.0));
  Spack result;
  const auto tmelt = C::Tmelt;
  Smask ice_mask = (t < tmelt) && ice;
  Smask liq_mask = !ice_mask;

  // -------------------------------------------

  // Flatau formulation:
  if (ice_mask.any()) {
    Spack ice_result = (ai[0] + dt*(ai[1]+dt*(ai[2]+dt*(ai[3]+dt*(ai[4]+dt*(ai[5]+dt*(ai[6]+dt*(ai[7]+
                                                                                                ai[8]*dt))))))))*100;
    result.set(ice_mask, ice_result);
  }
  if (liq_mask.any()) {
    Spack liq_result = (a[0] + dt*(a[1]+dt*(a[2]+dt*(a[3]+dt*(a[4]+dt*(a[5]+dt*(a[6]+dt*(a[7]+a[8]*dt))))))))*100;
    result.set(liq_mask, liq_result);
  }

  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::qv_sat(const Spack& t_atm, const Spack& p_atm, const bool ice)
{
  Spack e_pres; // saturation vapor pressure [Pa]

  //e_pres = polysvp1(t_atm, ice);
  e_pres = svp_murphy_koop(t_atm, ice);
  const auto ep_2 = C::ep_2;
  return ep_2 * e_pres / pack::max(p_atm-e_pres, sp(1.e-3));
}

} // namespace physics
} // namespace scream

#endif
