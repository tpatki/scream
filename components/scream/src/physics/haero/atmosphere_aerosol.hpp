#ifndef SCREAM_HIRES_AEROSOL_HPP
#define SCREAM_HIRES_AEROSOL_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
/*
 * The class responsible to handle the atmospheric aerosols
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
*/
{
class HiResAerosol : public AtmosphereProcess
{
public:
  using field_type       = Field<Real,       device_type>;
  using const_field_type = Field<const Real, device_type>;

  // Constructors
  HiResAerosol(const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const {return AtmosphereProcessType::Physics;}

  // The name of the subcomponent
  std::string name () const { return "Aerosol Microphysics"; }

  // The communicator used by subcomponent
  const ekat::Comm& get_comm () const { return m_p3_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_p3_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // The three main interfaces for the subcomponent
  void initialize (const util::TimeStamp& t0);
  void run        (const Real dt);
  void finalize   ();

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const;

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>& get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>& get_computed_fields () const { return m_computed_fields; }

protected:
  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real, device_type>& f);
  void set_computed_field_impl (const Field<      Real, device_type>& f);

  std::map<std::string,const_field_type>  m_haero_fields_in;
  std::map<std::string,field_type>        m_haero_fields_out;

  std::map<std::string,const_field_type::view_type::HostMirror>  m_haero_host_views_in;
  std::map<std::string,field_type::view_type::HostMirror>        m_haero_host_views_out;

  // Used to init some fields. For now, only needed for stand-alone haero runs
  std::shared_ptr<FieldInitializer>  m_initializer;

  util::TimeStamp     m_current_ts;
  ekat::Comm          m_haero_comm;

  ekat::ParameterList m_haero_params;
};

}
