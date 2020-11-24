#ifndef SCREAM_P3_INTERFACE_HPP
#define SCREAM_P3_INTERFACE_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_assert.hpp"

// Put everything into a scream namespace
namespace scream {

extern "C"
{

// Fortran routines to be called from C
void p3_init_f90 (Int num_cols, Int num_levs);
void p3_standalone_init_f90 (Real* q, Real* T_atm, Real* zi, Real* pmid, Real* dpres,
                             Real* ast, Real* ni_activated, Real* nc_nuceat_tend, Real* qv_prev, Real* T_prev,
                             Real* qv, Real* qc, Real* qr, Real* qi, Real* qm, Real* nc, Real* nr, Real* ni, Real* bm
                            );
void p3_main_f90 (const Real& dtime,
                  const Real* zi, const Real* pmid,
                  const Real* dpres, const Real* ast,
                  const Real* ni_activated, const Real* nc_nuceat_tend,
                  Real* q, Real* FQ, Real* T, Real* qv_prev, Real* T_prev,
                  const Real* exner, const Real* dz);
void p3_finalize_f90 ();

void p3_main_c2f(
         Real* qc,
         Real* nc,
         Real* qr,
         Real* nr,
         Real* th_atm,
         Real* qv,
         const Real& dtime,
         Real* qi,
         Real* qm,
         Real* ni,
         Real* bm,
         Real* pmid              ,
         Real* dz                ,
         Real* nc_nuceat_tend    ,
         Real* nccn_prescribed   ,
         Real* ni_activated      ,
         Real* inv_qc_relvar     ,
         const Int& it,                
         Real* precip_liq_surf   ,
         Real* precip_ice_surf   ,
         const Int& its,               
         const Int& ite,               
         const Int& kts,               
         const Int& kte,               
         Real* rel               ,
         Real* rei               ,
         Real* rho_qi            ,
         Real* dpres             ,
         Real* exner             ,
         Real* qv2qi_depos_tend  ,
         Real* precip_total_tend ,
         Real* nevapr            ,
         Real* qr_evap_tend      ,
         Real* precip_liq_flux   ,
         Real* precip_ice_flux   ,
         Real* cld_frac_r        ,
         Real* cld_frac_l        ,
         Real* cld_frac_i        ,
         Real* mu                ,
         Real* lambdac           ,
         Real* liq_ice_exchange  ,
         Real* vap_liq_exchange  ,
         Real* vap_ice_exchange  ,
         Real* qv_prev           ,
         Real* t_prev            ,
         Real& elapsed_s 
         );

} // extern "C"

} // namespace scream

#endif // SCREAM_P3_INTERFACE_HPP
