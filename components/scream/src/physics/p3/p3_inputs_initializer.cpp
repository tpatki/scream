#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/p3/scream_p3_interface.hpp"

#include <array>

namespace scream
{

void P3InputsInitializer::add_field (const field_type &f)
{
  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}

void P3InputsInitializer::
add_field (const field_type &f, const field_type& f_ref,
           const remapper_ptr_type& remapper)
{
  if (m_remapper) {
    // Sanity check
    EKAT_REQUIRE_MSG (m_remapper->get_src_grid()->name()==remapper->get_src_grid()->name(),
      "Error! A remapper was already set in P3InputsInitializer, but its src grid differs from"
      "       the grid of the input remapper of this call.\n");
  } else {
    m_remapper = remapper;
    m_remapper->registration_begins();
  }

  const auto& id = f.get_header().get_identifier();
  const auto& id_ref = f_ref.get_header().get_identifier();

  // To the AD, we only expose the fact that we init f_ref...
  m_fields_id.insert(id_ref);

  // ...but P3 only knows how to init f...
  m_fields.emplace(id.name(),f);

  // ...hence, we remap to f_ref.
  m_remapper->register_field(f, f_ref);
}


// =========================================================================================
void P3InputsInitializer::initialize_fields ()
{
  // Safety check: if we're asked to init anything at all,
  // then we should have been asked to init 10 fields.
  int count = 0;
  count += m_fields.count("q");
  count += m_fields.count("T");
  count += m_fields.count("ast");
  count += m_fields.count("ni_activated");
  count += m_fields.count("nc_nuceat_tend");
  count += m_fields.count("pmid");
  count += m_fields.count("dp");
  count += m_fields.count("zi");
  count += m_fields.count("qv_prev");
  count += m_fields.count("T_prev");
  
  count += m_fields.count("qv");
  count += m_fields.count("qc");
  count += m_fields.count("qr");
  count += m_fields.count("qi");
  count += m_fields.count("qm");
  count += m_fields.count("nc");
  count += m_fields.count("nr");
  count += m_fields.count("ni");
  count += m_fields.count("bm");
  if (count==0) {
    return;
  }

  EKAT_REQUIRE_MSG (count==19,
    "Error! P3InputsInitializer is expected to init 'q','T','ast','ni_activated','nc_nuceat_tend','pmid','dp','zi','qv_prev','T_prev',.\n"
    "       'qv', 'qc', 'qr', 'qi', 'qm', 'nc', 'nr', 'ni', 'bm'.\n"
    "       Only " + std::to_string(count) + " of those have been found.\n"
    "       Please, check the atmosphere processes you are using,"
    "       and make sure they agree on who's initializing each field.\n");

  // Get device views
  auto d_q     = m_fields.at("q").get_view();
  auto d_T     = m_fields.at("T").get_view();
  auto d_ast   = m_fields.at("ast").get_view();
  auto d_ni_activated  = m_fields.at("ni_activated").get_view();
  auto d_nc_nuceat_tend = m_fields.at("nc_nuceat_tend").get_view();
  auto d_pmid  = m_fields.at("pmid").get_view();
  auto d_dpres  = m_fields.at("dp").get_view();
  auto d_zi    = m_fields.at("zi").get_view();
  auto d_qv_prev = m_fields.at("qv_prev").get_view();
  auto d_t_prev  = m_fields.at("T_prev").get_view();
  
  auto d_qv    = m_fields.at("qv").get_view();
  auto d_qc    = m_fields.at("qc").get_view();
  auto d_qr    = m_fields.at("qr").get_view();
  auto d_qi    = m_fields.at("qi").get_view();
  auto d_qm    = m_fields.at("qm").get_view();
  auto d_nc    = m_fields.at("nc").get_view();
  auto d_nr    = m_fields.at("nr").get_view();
  auto d_ni    = m_fields.at("ni").get_view();
  auto d_bm    = m_fields.at("bm").get_view();
  // Create host mirrors
  auto h_q     = Kokkos::create_mirror_view(d_q);
  auto h_T     = Kokkos::create_mirror_view(d_T);
  auto h_ast   = Kokkos::create_mirror_view(d_ast);
  auto h_ni_activated  = Kokkos::create_mirror_view(d_ni_activated);
  auto h_nc_nuceat_tend = Kokkos::create_mirror_view(d_nc_nuceat_tend);
  auto h_pmid  = Kokkos::create_mirror_view(d_pmid);
  auto h_dpres  = Kokkos::create_mirror_view(d_dpres);
  auto h_zi    = Kokkos::create_mirror_view(d_zi);
  auto h_qv_prev = Kokkos::create_mirror_view(d_qv_prev);
  auto h_t_prev = Kokkos::create_mirror_view(d_t_prev);
  
  auto h_qv     = Kokkos::create_mirror_view(d_qv);
  auto h_qc     = Kokkos::create_mirror_view(d_qc);
  auto h_qr     = Kokkos::create_mirror_view(d_qr);
  auto h_qi     = Kokkos::create_mirror_view(d_qi);
  auto h_qm     = Kokkos::create_mirror_view(d_qm);
  auto h_nc     = Kokkos::create_mirror_view(d_nc);
  auto h_nr     = Kokkos::create_mirror_view(d_nr);
  auto h_ni     = Kokkos::create_mirror_view(d_ni);
  auto h_bm     = Kokkos::create_mirror_view(d_bm);
  // Get host mirros' raw pointers
  auto q     = h_q.data();
  auto T_atm     = h_T.data();
  auto ast   = h_ast.data();
  auto ni_activated  = h_ni_activated.data();
  auto nc_nuceat_tend = h_nc_nuceat_tend.data();
  auto pmid  = h_pmid.data();
  auto dpres  = h_dpres.data();
  auto zi    = h_zi.data();
  auto qv_prev = h_qv_prev.data();
  auto t_prev  = h_t_prev.data();
  
  auto qv = h_qv.data();
  auto qc = h_qc.data();
  auto qr = h_qr.data();
  auto qi = h_qi.data();
  auto qm = h_qm.data();
  auto nc = h_nc.data();
  auto nr = h_nr.data();
  auto ni = h_ni.data();
  auto bm = h_bm.data();
  // Call f90 routine
  p3_standalone_init_f90 (q, T_atm, zi, pmid, dpres, ast, ni_activated, nc_nuceat_tend, qv_prev, t_prev,
                          qv, qc, qr, qi, qm, nc, nr, ni, bm);

  // Deep copy back to device
  Kokkos::deep_copy(d_q,h_q);
  Kokkos::deep_copy(d_T,h_T);
  Kokkos::deep_copy(d_ast,h_ast);
  Kokkos::deep_copy(d_ni_activated,h_ni_activated);
  Kokkos::deep_copy(d_nc_nuceat_tend,h_nc_nuceat_tend);
  Kokkos::deep_copy(d_pmid,h_pmid);
  Kokkos::deep_copy(d_dpres,h_dpres);
  Kokkos::deep_copy(d_zi,h_zi);
  Kokkos::deep_copy(d_qv_prev,h_qv_prev);
  Kokkos::deep_copy(d_t_prev,h_t_prev);

  Kokkos::deep_copy(d_qv,h_qv);
  Kokkos::deep_copy(d_qc,h_qc);
  Kokkos::deep_copy(d_qr,h_qr);
  Kokkos::deep_copy(d_qi,h_qi);
  Kokkos::deep_copy(d_qm,h_qm);
  Kokkos::deep_copy(d_nc,h_nc);
  Kokkos::deep_copy(d_nr,h_nr);
  Kokkos::deep_copy(d_ni,h_ni);
  Kokkos::deep_copy(d_bm,h_bm);
  // If we are in charge of init-ing FQ as well, init it to 0.
  if (m_fields.count("FQ")==1) {
    // Init FQ to 0
    auto d_FQ = m_fields.at("FQ").get_view();
    Kokkos::deep_copy(d_FQ,Real(0));
  }

  if (m_remapper) {
    m_remapper->registration_ends();

    m_remapper->remap(true);

    // Now we can destroy the remapper
    m_remapper = nullptr;
  }
}

} // namespace scream
