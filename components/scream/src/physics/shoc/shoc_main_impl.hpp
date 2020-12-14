#ifndef SHOC_MAIN_IMPL_HPP
#define SHOC_MAIN_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_main. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_main(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Int&                   num_qtracers,
  const Int&                   nadv,
  const Scalar&                dtime,
  const Scalar&                host_dx,
  const Scalar&                host_dy,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& presi,
  const uview_1d<const Spack>& pdel,
  const uview_1d<const Spack>& thv,
  const uview_1d<const Spack>& w_field,
  const Scalar&                wthl_sfc,
  const Scalar&                wqw_sfc,
  const Scalar&                uw_sfc,
  const Scalar&                vw_sfc,
  const uview_1d<const Spack>& wtracer_sfc,
  const uview_1d<const Spack>& exner,
  const Scalar&                phis,
  const uview_1d<Spack>&       rho_zt,
  const uview_1d<Spack>&       shoc_qv,
  const uview_1d<Spack>&       dz_zt,
  const uview_1d<Spack>&       dz_zi,

  const uview_1d<Spack>&       thv_zi,
  const uview_1d<Spack>&       sterm,
  const uview_1d<Spack>&       sterm_zt,
  const uview_1d<Spack>&       a_diss,
  const uview_1d<Spack>&       rdp_zt,
  const uview_1d<Spack>&       tmpi,
  const uview_1d<Spack>&       tkh_zi,
  const uview_1d<Spack>&       tk_zi,
  const uview_1d<Spack>&       rho_zi,
  const uview_1d<Scalar>&      du,
  const uview_1d<Scalar>&      dl,
  const uview_1d<Scalar>&      d,
  const uview_2d<Spack>&       X1,
  const uview_1d<Spack>&       isotropy_zi,
  const uview_1d<Spack>&       w_sec_zi,
  const uview_1d<Spack>&       brunt_zi,
  const uview_1d<Spack>&       thetal_zi,
  const uview_1d<Spack>&       wthl_sec_zt,
  const uview_1d<Spack>&       wqw_sec_zt,
  const uview_1d<Spack>&       w3_zt,
  const uview_1d<Spack>&       thl_sec_zt,
  const uview_1d<Spack>&       qwthl_sec_zt,
  const uview_1d<Spack>&       qw_sec_zt,

  const uview_1d<Spack>&       host_dse,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       thetal,
  const uview_1d<Spack>&       qw,
  const uview_1d<Spack>&       u_wind,
  const uview_1d<Spack>&       v_wind,
  const uview_1d<Spack>&       wthv_sec,
  const uview_2d<Spack>&       qtracers,
  const uview_1d<Spack>&       tk,
  const uview_1d<Spack>&       tkh,
  const uview_1d<Spack>&       shoc_cldfrac,
  const uview_1d<Spack>&       shoc_ql,
  Scalar&                      pblh,
  const uview_1d<Spack>&       shoc_ql2,
  const uview_1d<Spack>&       shoc_mix,
  const uview_1d<Spack>&       w_sec,
  const uview_1d<Spack>&       thl_sec,
  const uview_1d<Spack>&       qw_sec,
  const uview_1d<Spack>&       qwthl_sec,
  const uview_1d<Spack>&       wthl_sec,
  const uview_1d<Spack>&       wqw_sec,
  const uview_1d<Spack>&       wtke_sec,
  const uview_1d<Spack>&       uw_sec,
  const uview_1d<Spack>&       vw_sec,
  const uview_1d<Spack>&       w3,
  const uview_1d<Spack>&       wqls_sec,
  const uview_1d<Spack>&       brunt,
  const uview_1d<Spack>&       isotropy)
{
  Scalar ustar(0), kbfs(0), obklen(0);

  // Compute integrals of static energy, kinetic energy, water vapor, and liquid water
  // for the computation of total energy before SHOC is called.  This is for an
  // effort to conserve energy since liquid water potential temperature (which SHOC
  // conserves) and static energy (which E3SM conserves) are not exactly equal.
  Scalar se_b(0), ke_b(0), wv_b(0), wl_b(0);
  shoc_energy_integrals(team,nlev,host_dse,pdel,qw,shoc_ql,u_wind,v_wind,
                        se_b,ke_b,wv_b,wl_b);

  for (Int t=0; t<nadv; ++t) {
    // Check TKE to make sure values lie within acceptable
    // bounds after host model performs horizontal advection
    check_tke(team,nlev,tke);

    // Define vertical grid arrays needed for
    // vertical derivatives in SHOC, also
    // define air density
    shoc_grid(team,nlev,nlevi,zt_grid,zi_grid,pdel,dz_zt,dz_zi,rho_zt);

    // Compute the planetary boundary layer height, which is an
    // input needed for the length scale calculation.

    // Update SHOC water vapor, to be used by the next two routines
    compute_shoc_vapor(team,nlev,qw,shoc_ql,shoc_qv);

    const Int nlev_view_indx = nlev/Spack::n;
    const Int nlev_pack_indx = nlev%Spack::n;
    shoc_diag_obklen(uw_sfc,vw_sfc,wthl_sfc,
                     wqw_sfc,thetal(nlev_view_indx)[nlev_pack_indx],
                     shoc_ql(nlev_view_indx)[nlev_pack_indx],
                     shoc_qv(nlev_view_indx)[nlev_pack_indx],
                     ustar,kbfs,obklen);


//    call pblintd(&
//       shcol,nlev,nlevi,&                   ! Input
//       zt_grid,zi_grid,thetal,shoc_ql,&     ! Input
//       shoc_qv,u_wind,v_wind,&              ! Input
//       ustar,obklen,kbfs,shoc_cldfrac,&     ! Input
//       pblh)                                ! Output

    // Update the turbulent length scale
    shoc_length(team,nlev,nlevi,host_dx,host_dy,
                pblh,tke,zt_grid,zi_grid,dz_zt,
                dz_zi,wthv_sec,thetal,thv,thv_zi,
                brunt,shoc_mix);

    // Advance the SGS TKE equation
    shoc_tke(team,nlev,nlevi,dtime,wthv_sec,shoc_mix,
             dz_zi,dz_zt,pres,u_wind,v_wind,brunt,obklen,
             zt_grid,zi_grid,pblh,sterm,sterm_zt,a_diss,tke,tk,tkh,isotropy);

    // Update SHOC prognostic variables here
    // via implicit diffusion solver
    update_prognostics_implicit(team,nlev,nlevi,num_qtracers,
                                dtime,dz_zt,dz_zi,rho_zt,
                                zt_grid,zi_grid,tk,tkh,
                                uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,
                                wtracer_sfc,
                                rdp_zt,tmpi,tkh_zi,tk_zi,rho_zi,du,dl,d,X1,
                                thetal,
                                qw,qtracers,tke,u_wind,v_wind);

    // Diagnose the second order moments
    Scalar ustar2(0), wstar(0);
    diag_second_shoc_moments(team,nlev,nlevi,thetal,qw,u_wind,v_wind,
                             tke,isotropy,tkh,tk,dz_zi,zt_grid,zi_grid,
                             shoc_mix,wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,
                             ustar2,wstar,isotropy_zi,tkh_zi,tk_zi,
                             thl_sec,qw_sec,wthl_sec,wqw_sec,qwthl_sec,
                             uw_sec,vw_sec,wtke_sec,w_sec);

    // Diagnose the third moment of vertical velocity,
    //  needed for the PDF closure
    diag_third_shoc_moments(team,nlev,nlevi,w_sec,thl_sec,wthl_sec,
                            isotropy,brunt,thetal,tke,dz_zt,dz_zi,
                            zt_grid,zi_grid,w_sec_zi,isotropy_zi,
                            brunt_zi,thetal_zi,w3);

    // Call the PDF to close on SGS cloud and turbulence
    shoc_assumed_pdf(team,nlev,nlevi,thetal,qw,w_field,thl_sec,
                     qw_sec,wthl_sec,w_sec,wqw_sec,qwthl_sec,w3,
                     pres,zt_grid,zi_grid,

                     wthl_sec_zt,wqw_sec_zt,w3_zt,thl_sec_zt,qwthl_sec_zt,qw_sec_zt,

                     shoc_cldfrac,shoc_ql,wqls_sec,wthv_sec,shoc_ql2);

    // Check TKE to make sure values lie within acceptable
    // bounds after vertical advection, etc.
    check_tke(team,nlev,tke);
  }

  // End SHOC parameterization

  // Use SHOC outputs to update the host model
  // temperature
  update_host_dse(team,nlev,thetal,shoc_ql,exner,zt_grid,
                  phis,host_dse);

  Scalar se_a(0), ke_a(0), wv_a(0), wl_a(0);
  shoc_energy_integrals(team,nlev,host_dse,pdel,qw,shoc_ql,u_wind,v_wind,
                        se_a,ke_a,wv_a,wl_a);

  shoc_energy_fixer(team,nlev,nlevi,dtime,nadv,zt_grid,zi_grid,
                    se_b,ke_b,wv_b,wl_b,se_a,ke_a,wv_a,wl_a,
                    wthl_sfc,wqw_sfc,rho_zt,tke,presi,rho_zi,host_dse);

  // Remaining code is to diagnose certain quantities
  // related to PBL.  No answer changing subroutines
  // should be placed at this point onward.

  // Update PBLH, as other routines outside of SHOC
  // may require this variable.

  // Update SHOC water vapor, to be used by the next two routines
  compute_shoc_vapor(team,nlev,qw,shoc_ql,shoc_qv);

  const Int nlev_view_indx = nlev/Spack::n;
  const Int nlev_pack_indx = nlev%Spack::n;
  shoc_diag_obklen(uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,
                   thetal(nlev_view_indx)[nlev_pack_indx],
                   shoc_ql(nlev_view_indx)[nlev_pack_indx],
                   shoc_qv(nlev_view_indx)[nlev_pack_indx],
                   ustar,kbfs,obklen);

//  pblintd(&
//     shcol,nlev,nlevi,&                   ! Input
//     zt_grid,zi_grid,thetal,shoc_ql,&     ! Input
//     shoc_qv,u_wind,v_wind,&              ! Input
//     ustar,obklen,kbfs,shoc_cldfrac,&     ! Input
//     pblh)                                ! Output
  pblh = 0;

}

} // namespace shoc
} // namespace scream

#endif
