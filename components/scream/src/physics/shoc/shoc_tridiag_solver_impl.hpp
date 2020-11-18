#ifndef SHOC_TRIDIAG_SOLVER_IMPL_HPP
#define SHOC_TRIDIAG_SOLVER_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "ekat/util/ekat_tridiag.hpp"

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_shoc_decomp(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& kv_term,
  const uview_1d<const Spack>& tmpi,
  const uview_1d<const Spack>& rdp_zt,
  const Scalar&                dtime,
  const Scalar&                flux,
  const uview_1d<Spack>&       du,
  const uview_1d<Spack>&       dl,
  const uview_1d<Spack>&       d)
{
  const auto ggr = C::gravit;

  const auto skv_term = scalarize(kv_term);
  const auto stmpi = scalarize(tmpi);

  const Int nlev_pack = ekat::npack<Spack>(nlev);

  // Compute entries of the tridiagonal system
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {
    // Compute shift of kv_term and tmpi
    Spack kv_term_k, kv_term_kp1, tmpi_k, tmpi_kp1;
    auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    auto shift_range = range_pack;
    shift_range.set(range_pack > nlev-1, 1); // don't calculate shift above nlev-1
    ekat::index_and_shift<1>(skv_term, shift_range, kv_term_k, kv_term_kp1);
    ekat::index_and_shift<1>(stmpi, shift_range, tmpi_k, tmpi_kp1);

    // Determine superdiagonal (du) and subdiagonal (dl) coeffs of the
    // tridiagonal diffusion matrix.
    du(k) = -kv_term_kp1*tmpi_kp1*rdp_zt(k);
    dl(k) = -kv_term(k)*tmpi(k)*rdp_zt(k);

    // The bottom element of the superdiagonal (du) and the top element of
    // the subdiagonal (dl) is set to zero (not included in linear system).
    du(k).set(range_pack == nlev-1, 0);
    dl(k).set(range_pack == 0, 0);

    // The diagonal elements are a combination of du and dl (d=1-du-dl). Surface
    // fluxes are applied explicitly in the diagonal at the top level.
    d(k).set(range_pack == 0, 1 - du(k));
    d(k).set(0 < range_pack && range_pack < nlev-1, 1 - du(k) - dl(k));
    d(k).set(range_pack == nlev-1, 1 - dl(k) + flux*dtime*ggr*rdp_zt(k));
  });
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_shoc_solve(
  const MemberType&      team,
  const Int&             nlev,
  const Int&             num_rhs,
  const uview_1d<const Spack>& du,
  const uview_1d<const Spack>& dl,
  const uview_1d<const Spack>& d,
  const uview_2d<Spack>&       var)
{
//  for (unsigned int p=0; p<num_rhs; ++p) {
//  for(unsigned int k=0; k<nlev; ++k) {
//      std::cout << team.league_rank()+1 << "," << k+1 << "," << p+1 << ":   "
//                << var(k,p/Spack::n)[p%Spack::n] << std::endl;
//    }
//  std::cout << std::endl;
//  }

  Kokkos::View<Scalar*, Kokkos::LayoutRight>
    du_in("du", nlev),
    dl_in("dl", nlev),
    d_in ("d",  nlev);
  for (Int k=0; k<nlev; ++k) {
    const auto view_indx = k/Spack::n;
    const auto pack_indx = k%Spack::n;

    du_in(k) = du(view_indx)[pack_indx];
    dl_in(k) = dl(view_indx)[pack_indx];
    d_in(k)  = d (view_indx)[pack_indx];
  }



  ekat::bfb(team,dl_in,d_in,du_in,var);

}

} // namespace shoc
} // namespace scream

#endif
