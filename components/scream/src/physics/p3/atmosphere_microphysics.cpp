#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/p3/p3_f90.hpp"
#include "physics/p3/p3_main_impl.hpp"

#include "physics/share/physics_constants.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <array>

namespace scream
{
/*
 * P3 Microphysics routines
*/

// =========================================================================================
P3Microphysics::P3Microphysics (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_p3_comm (comm)
 , m_p3_params (params)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, table lookups, p3 options.
*/
  m_initializer = create_field_initializer<P3InputsInitializer>();
}

// =========================================================================================
void P3Microphysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");

  constexpr int NVL = SCREAM_NUM_VERTICAL_LEV;
  constexpr int QSZ =  35;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

  const auto& grid_name = m_p3_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  const int num_dofs = grid->get_num_local_dofs();
  const int nc = num_dofs;

  m_num_cols = nc;
  m_num_levs = NVL;

  using namespace ShortFieldTagsNames;

  FieldLayout scalar3d_layout_mid { {COL,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar2d_layout_mid { {COL}, {nc} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar3d_layout_int { {COL,VL}, {nc,NVL+1} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout_mid{ {COL,CMP,VL}, {nc,QSZ,NVL} };
  FieldLayout tracers_layout { {COL,VAR,VL}, {nc,QSZ,NVL} };

  // Inputs
  auto nondim = m/m;
  m_required_fields.emplace("ast",            scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("pmid",           scalar3d_layout_mid,   Pa, grid_name);
//  m_required_fields.emplace("dp",             scalar3d_layout_mid,   Pa, grid_name);
  m_required_fields.emplace("zi",             scalar3d_layout_int,   m, grid_name);
//  m_required_fields.emplace("qv_prev", vector3d_layout_mid, Q, grid_name);
  // Input-Outputs
  m_required_fields.emplace("FQ", tracers_layout,      Q, grid_name);
  m_required_fields.emplace("T",  scalar3d_layout_mid, K, grid_name);
  m_required_fields.emplace("T_atm",  scalar3d_layout_mid, K, grid_name);
  m_required_fields.emplace("q",  vector3d_layout_mid, Q, grid_name);

  m_computed_fields.emplace("FQ", tracers_layout,      Q, grid_name);
  m_computed_fields.emplace("T",  scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("T_atm",  scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("pmid",           scalar3d_layout_mid,   Pa, grid_name);
  m_computed_fields.emplace("zi",             scalar3d_layout_int,   m, grid_name);
  m_computed_fields.emplace("ast",            scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("q",  vector3d_layout_mid, Q, grid_name);
//  m_computed_fields.emplace("qv_prev", vector3d_layout_mid, Q, grid_name);

  // TODO: Fix the units for these variables
  // Prognostic State
  m_required_fields.emplace("qv",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("qc",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("qr",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("qi",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("qm",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("nc",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("nr",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("ni",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("bm",             scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("th_atm",         scalar3d_layout_mid,   nondim, grid_name);
  //
  m_computed_fields.emplace("qv",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("qc",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("qr",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("qi",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("qm",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("nc",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("nr",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("ni",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("bm",             scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("th_atm",         scalar3d_layout_mid,   nondim, grid_name);
  // Diagnostic Inputs
  m_required_fields.emplace("nc_nuceat_tend",  scalar3d_layout_mid,   1/(kg*s), grid_name);
  m_required_fields.emplace("nccn_prescribed", scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("ni_activated",    scalar3d_layout_mid,   1/kg, grid_name);
  m_required_fields.emplace("inv_qc_relvar",   scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("cld_frac_i",      scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("cld_frac_l",      scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("cld_frac_r",      scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("pres",            scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("dz",              scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("dp",              scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("exner",           scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("qv_prev",         scalar3d_layout_mid,        Q, grid_name);
  m_required_fields.emplace("T_prev",          scalar3d_layout_mid,        K, grid_name);
  //
  m_computed_fields.emplace("nc_nuceat_tend",  scalar3d_layout_mid,   1/(kg*s), grid_name);
  m_computed_fields.emplace("nccn_prescribed", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("ni_activated",    scalar3d_layout_mid,   1/kg, grid_name);
  m_computed_fields.emplace("inv_qc_relvar",   scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("cld_frac_i",      scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("cld_frac_l",      scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("cld_frac_r",      scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("pres",            scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("dz",              scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("dp",              scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("exner",           scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("qv_prev",         scalar3d_layout_mid,        Q, grid_name);
  m_computed_fields.emplace("T_prev",          scalar3d_layout_mid,        K, grid_name);
  // Diagnostic Outputs
  m_computed_fields.emplace("mu_c",               scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("lamc",               scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("qv2qi_depos_tend",   scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("precip_liq_flux",    scalar3d_layout_int,   nondim, grid_name);
  m_computed_fields.emplace("precip_ice_flux",    scalar3d_layout_int,   nondim, grid_name);
  m_computed_fields.emplace("precip_liq_surf",    scalar2d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("precip_ice_surf",    scalar2d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("diag_eff_radius_qc", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("diag_eff_radius_qi", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("rho_qi", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("precip_total_tend", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("nevapr", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("qr_evap_tend", scalar3d_layout_mid,   nondim, grid_name);
  // History Only
  m_computed_fields.emplace("liq_ice_exchange", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("vap_liq_exchange", scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("vap_ice_exchange", scalar3d_layout_mid,   nondim, grid_name);
}

// =========================================================================================
void P3Microphysics::initialize_impl (const util::TimeStamp& t0)
{
  using namespace p3;
  m_current_ts = t0;

  // Call f90 routine
//  p3_init_f90 (m_num_cols,m_num_levs);
  p3_init();

  // We may have to init some fields from within P3. This can be the case in a P3 standalone run.
  // Some options:
  //  - we can tell P3 it can init all inputs or specify which ones it can init. We call the
  //    resulting list of inputs the 'initializaable' (or initable) inputs. The default is
  //    that no inputs can be inited.
  //  - we can request that P3 either inits no inputs or all of the initable ones (as specified
  //    at the previous point). The default is that P3 must be in charge of init ing ALL or NONE
  //    of its initable inputs.
  // Recall that:
  //  - initable fields may not need initialization (e.g., some other atm proc that
  //    appears earlier in the atm dag might provide them).

  std::vector<std::string> p3_inputs = {"q","T_atm","FQ","ast","ni_activated","nc_nuceat_tend","pmid","dp","zi","qv_prev","T_prev",
                                        "qv", "qc", "qr", "qi", "qm", "nc", "nr", "ni", "bm" 
                                       };
  using strvec = std::vector<std::string>;
  const strvec& allowed_to_init = m_p3_params.get<strvec>("Initializable Inputs",strvec(0));
  const bool can_init_all = m_p3_params.get<bool>("Can Initialize All Inputs", false);
  const bool init_all_or_none = m_p3_params.get<bool>("Must Init All Inputs Or None", true);

  const strvec& initable = can_init_all ? p3_inputs : allowed_to_init;
  if (initable.size()>0) {
    bool all_inited = true, all_uninited = true;
    for (const auto& name : initable) {
      const auto& f = m_p3_fields_in.at(name);
      const auto& track = f.get_header().get_tracking();
      if (track.get_init_type()==InitType::None) {
        // Nobody claimed to init this field. P3InputsInitializer will take care of it
        m_initializer->add_me_as_initializer(f);
        all_uninited &= true;
        all_inited &= false;
      } else {
        all_uninited &= false;
        all_inited &= true;
      }
    }

    // In order to gurantee some consistency between inputs, it is best if P3
    // initializes either none or all of the inputs.
    EKAT_REQUIRE_MSG (!init_all_or_none || all_inited || all_uninited,
                      "Error! Some p3 inputs were marked to be inited by P3, while others weren't.\n"
                      "       P3 was requested to init either all or none of the inputs.\n");
  }
}

// =========================================================================================
void P3Microphysics::run_impl (const Real dt)
{
  using namespace p3;
  // std::array<const char*, num_views> view_names = {"q", "FQ", "T", "zi", "pmid", "dpres", "ast", "ni_activated", "nc_nuceat_tend"};

  std::vector<const Real*> in;
  std::vector<Real*> out;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_p3_fields_in) {
    Kokkos::deep_copy(m_p3_host_views_in.at(it.first),it.second.get_view());
  }
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(m_p3_host_views_out.at(it.first),it.second.get_view());
  }

  // Deal with local arrays that define inputs to p3_main
  // TODO: Some of these views should really just be input.  but at the moment I can't figure out
  // how to pass an input only view to a subroutine without getting a build error.  So I've made
  // every field an input/output for the time being until I can cross this bridge - Aaron
  using P3F  = Functions<Real, DefaultDevice>;
  using PC = scream::physics::Constants<Real>;
//  using namespace ekat::Pack;
  using Spack  = typename P3F::Spack;
  using Pack   = typename ekat::Pack<Real, 1>;

  auto pmid   = m_p3_fields_out["pmid"].get_reshaped_view<Pack**>();
  auto T_atm  = m_p3_fields_out["T_atm"].get_reshaped_view<Pack**>();
  auto zi     = m_p3_fields_out["zi"].get_reshaped_view<Pack**>();
  auto exner  = m_p3_fields_out["exner"].get_reshaped_view<Pack**>();
  auto dz     = m_p3_fields_out["dz"].get_reshaped_view<Pack**>();
  auto th_atm = m_p3_fields_out["th_atm"].get_reshaped_view<Pack**>();
  auto mu_c   = m_p3_fields_out["mu_c"].get_reshaped_view<Pack**>();
  auto lamc   = m_p3_fields_out["lamc"].get_reshaped_view<Pack**>();
  auto cld_frac_i = m_p3_fields_out["cld_frac_i"].get_reshaped_view<Pack**>();
  auto cld_frac_r = m_p3_fields_out["cld_frac_r"].get_reshaped_view<Pack**>();
  auto cld_frac_l = m_p3_fields_out["cld_frac_l"].get_reshaped_view<Pack**>();
  auto ast   = m_p3_fields_out["ast"].get_reshaped_view<Pack**>();
  auto qr    = m_p3_fields_out["qr"].get_reshaped_view<Pack**>();
  auto qi    = m_p3_fields_out["qi"].get_reshaped_view<Pack**>();
  auto inv_qc_relvar = m_p3_fields_out["inv_qc_relvar"].get_reshaped_view<Pack**>();
  Real mucon = 5.3;
  Real dcon  = 25.0 * pow(10.0,-6);
  Real qsmall = pow(10.0,-14); 
  Real mincld = 0.0001;  // TODO: These should be stored somewhere as more universal constants.  Or maybe in the P3 class hpp
  for (int i_col=0;i_col<m_num_cols;i_col++)
  {
    for (int i_lev=0;i_lev<m_num_levs;i_lev++)
    {
      exner(i_col,i_lev)  = 1.0/( pow( pmid(i_col,i_lev)*pow(10.0,-5), PC::RD*PC::INV_CP ) );
      th_atm(i_col,i_lev) = T_atm(i_col,i_lev) * exner(i_col,i_lev);
      dz(i_col,i_lev)     = zi(i_col,i_lev)-zi(i_col,i_lev+1);
      mu_c(i_col,i_lev)   = mucon;
      lamc(i_col,i_lev)   = (mucon - 1.0)/dcon;
      inv_qc_relvar(i_col,i_lev) = 1.0;
      // cloud fraction - TODO, this should be made into a universal function
      cld_frac_i(i_col,i_lev) = 1.0; //std::min(ast(i_col,i_lev),mincld);
      cld_frac_l(i_col,i_lev) = 1.0; //std::min(ast(i_col,i_lev),mincld);
      cld_frac_r(i_col,i_lev) = 1.0; //std::min(ast(i_col,i_lev),mincld);
      // Hard-code as "max_overlap" for now.  TODO: make more general, this can certainly be done when this is made a unverisal function.
      //if (i_lev != 0)
      //{
      //  if (qr(i_col,i_lev-1)>=qsmall or qi(i_col,i_lev-1)>=qsmall)
      //  {
      //    cld_frac_r(i_col,i_lev) = 1.0; //std::max(ast(i_col,i_lev-1),cld_frac_r(i_col,i_lev));
      //  }
      //}
    }
  }
  auto qc    = m_p3_fields_out["qc"].get_reshaped_view<Pack**>();
  auto qm    = m_p3_fields_out["qm"].get_reshaped_view<Pack**>();
  auto qv    = m_p3_fields_out["qv"].get_reshaped_view<Pack**>();
  auto nr    = m_p3_fields_out["nr"].get_reshaped_view<Pack**>();
  auto ni    = m_p3_fields_out["ni"].get_reshaped_view<Pack**>();
  auto nc    = m_p3_fields_out["nc"].get_reshaped_view<Pack**>();
  auto bm    = m_p3_fields_out["bm"].get_reshaped_view<Pack**>();

  auto nc_nuceat_tend  = m_p3_fields_out["nc_nuceat_tend"].get_reshaped_view<Pack**>();
  auto nccn_prescribed = m_p3_fields_out["nccn_prescribed"].get_reshaped_view<Pack**>();
  auto ni_activated    = m_p3_fields_out["ni_activated"].get_reshaped_view<Pack**>();
  auto dp              = m_p3_fields_out["dp"].get_reshaped_view<Pack**>();
  auto qv_prev         = m_p3_fields_out["qv_prev"].get_reshaped_view<Pack**>();
  auto T_prev          = m_p3_fields_out["T_prev"].get_reshaped_view<Pack**>();

  auto qv2qi_depos_tend   = m_p3_fields_out["qv2qi_depos_tend"].get_reshaped_view<Pack**>();
  auto diag_eff_radius_qc = m_p3_fields_out["diag_eff_radius_qc"].get_reshaped_view<Pack**>();
  auto diag_eff_radius_qi = m_p3_fields_out["diag_eff_radius_qi"].get_reshaped_view<Pack**>();
  auto rho_qi             = m_p3_fields_out["rho_qi"].get_reshaped_view<Pack**>();
  auto precip_total_tend  = m_p3_fields_out["precip_total_tend"].get_reshaped_view<Pack**>();
  auto nevapr             = m_p3_fields_out["nevapr"].get_reshaped_view<Pack**>();
  auto qr_evap_tend       = m_p3_fields_out["qr_evap_tend"].get_reshaped_view<Pack**>();
  auto precip_liq_surf    = m_p3_fields_out["precip_liq_surf"].get_view();
  auto precip_ice_surf    = m_p3_fields_out["precip_ice_surf"].get_view();
  auto precip_liq_flux    = m_p3_fields_out["precip_liq_flux"].get_reshaped_view<Pack**>();
  auto precip_ice_flux    = m_p3_fields_out["precip_ice_flux"].get_reshaped_view<Pack**>();

  using sview_2d   = typename KokkosTypes<DefaultDevice>::template view_2d<Real>;
  sview_2d col_location("col_location", m_num_cols, 3);


  auto liq_ice_exchange = m_p3_fields_out["liq_ice_exchange"].get_reshaped_view<Pack**>();
  auto vap_liq_exchange = m_p3_fields_out["vap_liq_exchange"].get_reshaped_view<Pack**>();
  auto vap_ice_exchange = m_p3_fields_out["vap_ice_exchange"].get_reshaped_view<Pack**>();
  // Pack our data into structs and ship it off to p3_main.
  m_it++;
  P3F::P3PrognosticState prog_state{
                                    qc,
                                    nc,
                                    qr,
                                    nr,
                                    qi,
                                    qm,
                                    ni,
                                    bm,
                                    qv,
                                    th_atm
                                   };
  P3F::P3DiagnosticInputs diag_inputs{
    nc_nuceat_tend,
    nccn_prescribed,
    ni_activated,
    inv_qc_relvar,
    cld_frac_i,
    cld_frac_l,
    cld_frac_r,
    pmid,
    dz,
    dp,
    exner,
    qv_prev,
    T_prev};
  P3F::P3DiagnosticOutputs diag_outputs{
    mu_c,
    lamc,
    qv2qi_depos_tend,
    precip_liq_surf,
    precip_ice_surf,
    diag_eff_radius_qc,
    diag_eff_radius_qi,
    rho_qi,
    precip_total_tend,
    nevapr,
    qr_evap_tend,
    precip_liq_flux,
    precip_ice_flux};
  P3F::P3Infrastructure infrastructure{dt, m_it, 1, m_num_cols, 1, m_num_levs,
                                       true, true, col_location};
  P3F::P3HistoryOnly history_only{
    liq_ice_exchange,
    vap_liq_exchange,
    vap_ice_exchange};

  using P3F  = Functions<Real, DefaultDevice>;
 
  auto q_before = qv(0,0);
  q_before = 0.0;
  for (int i_col=0;i_col<m_num_cols;i_col++)
  {
    for (int i_lev=0;i_lev<m_num_levs;i_lev++)
    {
      q_before = q_before + (qv(i_col,i_lev) + qc(i_col,i_lev) + qr(i_col,i_lev) + qi(i_col,i_lev) + qm(i_col,i_lev));
    }
  }
  auto elapsed_microsec = P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
                                       history_only, m_num_cols, m_num_levs);
  auto q_after = qv(0,0);
  q_after = 0.0; 
  for (int i_col=0;i_col<m_num_cols;i_col++)
  {
    for (int i_lev=0;i_lev<m_num_levs;i_lev++)
    {
      q_after = q_after + (qv(i_col,i_lev) + qc(i_col,i_lev) + qr(i_col,i_lev) + qi(i_col,i_lev) + qm(i_col,i_lev));
    }
  }
  printf("ASD = q_diff:  %f, %f, %f\n",q_before,q_after,q_after-q_before);
//  // Call f90 routine
//  Real elapsed_s;
//  p3_main_c2f(
//         m_raw_ptrs_out["qc"]                ,
//         m_raw_ptrs_out["nc"]                ,
//         m_raw_ptrs_out["qr"]                ,
//         m_raw_ptrs_out["nr"]                ,
//         m_raw_ptrs_out["th_atm"]            ,
//         m_raw_ptrs_out["qv"]                ,
//         dt,             
//         m_raw_ptrs_out["qi"]                ,
//         m_raw_ptrs_out["qm"]                ,
//         m_raw_ptrs_out["ni"]                ,
//         m_raw_ptrs_out["bm"]                ,
//         m_raw_ptrs_out["pmid"]              ,
//         m_raw_ptrs_out["dz"]                ,
//         m_raw_ptrs_out["nc_nuceat_tend"]    ,
//         m_raw_ptrs_out["nccn_prescribed"]   ,
//         m_raw_ptrs_out["ni_activated"]      ,
//         m_raw_ptrs_out["inv_qc_relvar"]     ,
//         m_it,                
//         m_raw_ptrs_out["precip_liq_surf"]   ,
//         m_raw_ptrs_out["precip_ice_surf"]   ,
//         1,               
//         m_num_cols,               
//         1,               
//         m_num_levs,               
//         m_raw_ptrs_out["diag_eff_radius_qc"]               ,
//         m_raw_ptrs_out["diag_eff_radius_qi"]               ,
//         m_raw_ptrs_out["rho_qi"]            ,
//         m_raw_ptrs_out["dp"]                ,
//         m_raw_ptrs_out["exner"]             ,
//         m_raw_ptrs_out["qv2qi_depos_tend"]  ,
//         m_raw_ptrs_out["precip_total_tend"] ,
//         m_raw_ptrs_out["nevapr"]            ,
//         m_raw_ptrs_out["qr_evap_tend"]      ,
//         m_raw_ptrs_out["precip_liq_flux"]   ,
//         m_raw_ptrs_out["precip_ice_flux"]   ,
//         m_raw_ptrs_out["cld_frac_r"]        ,
//         m_raw_ptrs_out["cld_frac_l"]        ,
//         m_raw_ptrs_out["cld_frac_i"]        ,
//         m_raw_ptrs_out["mu_c"]                ,
//         m_raw_ptrs_out["lamc"]           ,
//         m_raw_ptrs_out["liq_ice_exchange"]  ,
//         m_raw_ptrs_out["vap_liq_exchange"]  ,
//         m_raw_ptrs_out["vap_ice_exchange"]  ,
//         m_raw_ptrs_out["qv_prev"]           ,
//         m_raw_ptrs_out["T_prev"]            ,
//         elapsed_s); 
//  p3_main_f90 (dt, 
//               m_raw_ptrs_in["zi"], 
//               m_raw_ptrs_in["pmid"], 
//               m_raw_ptrs_in["dp"], 
//               m_raw_ptrs_in["ast"], 
//               m_raw_ptrs_in["ni_activated"], 
//               m_raw_ptrs_in["nc_nuceat_tend"], 
//               m_raw_ptrs_out["q"], 
//               m_raw_ptrs_out["FQ"], 
//               m_raw_ptrs_out["T_atm"], 
//               m_raw_ptrs_out["qv_prev"], 
//               m_raw_ptrs_out["T_prev"],
//               m_raw_ptrs_out["exner"],
//               m_raw_ptrs_out["dz"]
//              );

  // Copy outputs back to device
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(it.second.get_view(),m_p3_host_views_out.at(it.first));
  }

  // Get a copy of the current timestamp (at the beginning of the step) and
  // advance it, updating the p3 fields.
  auto ts = timestamp();
  ts += dt;
  m_p3_fields_out.at("q").get_header().get_tracking().update_time_stamp(ts);
  m_p3_fields_out.at("FQ").get_header().get_tracking().update_time_stamp(ts);
  m_p3_fields_out.at("T").get_header().get_tracking().update_time_stamp(ts);
}

// =========================================================================================
void P3Microphysics::finalize_impl()
{
//  p3_finalize_f90 ();
}

// =========================================================================================
void P3Microphysics::register_fields (FieldRepository<Real>& field_repo) const {
  for (auto& fid : m_required_fields) {
    field_repo.register_field(fid);
  }
  for (auto& fid : m_computed_fields) {
    field_repo.register_field(fid);
  }
}

void P3Microphysics::set_required_field_impl (const Field<const Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_in.emplace(name,f);
  m_p3_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_in[name] = m_p3_host_views_in[name].data();

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void P3Microphysics::set_computed_field_impl (const Field<      Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_out.emplace(name,f);
  m_p3_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_out[name] = m_p3_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
