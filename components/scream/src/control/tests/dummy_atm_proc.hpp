#include <catch2/catch.hpp>
#include "share/atmosphere_process.hpp"
#include "share/scream_pack.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/default_grid.hpp"
#include "control/atmosphere_driver.hpp"

namespace scream {

// === A dummy atm process, on Physics grid === //

template<typename DeviceType, int PackSize, bool forward>
class DummyProcess : public scream::AtmosphereProcess {
public:
  using device_type = DeviceType;

  DummyProcess (const Comm& comm, const ParameterList& params)
   : m_comm(comm)
  {
    m_params = params;
    m_id = comm.rank();

    if (forward) {
      m_name = "Physics_fwd";
    } else {
      m_name = "Physics_bwd";
    }
  }

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  std::set<std::string> get_required_grids () const {
    // TODO: define what grid the coupling runs on. Check with MOAB folks.
    static std::set<std::string> s;
    s.insert(m_name);
    return s;
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

  // The communicator associated with this atm process
  const Comm& get_comm () const { return m_comm; }

  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) {
    m_grid = grids_manager->get_grid(m_name);

    auto num_cols = m_grid->num_dofs();

    std::vector<FieldTag> tags = {FieldTag::Column,FieldTag::Component};
    std::vector<int> dims = {num_cols, m_params.get<int>("Number of vector components")};
    FieldLayout layout (tags,dims);

    std::string in_name = "field_";
    std::string out_name = "field_";
    if (forward) {
      in_name  += "0";
      out_name += "1";
    } else {
      in_name  += "1";
      out_name += "0";
    }

    m_input_fids.emplace(in_name,layout,units::one,m_grid->name());
    m_output_fids.emplace(out_name,layout,units::one,m_grid->name());
  }

  void initialize (const util::TimeStamp& t0) {
    m_time_stamp = t0;
  }

  void run (const double dt) {
    auto in = m_input.get_view();
    auto out = m_output.get_view();
    Kokkos::parallel_for(Kokkos::RangePolicy<>(0,in.size()),
      KOKKOS_LAMBDA(const int i) {
        out(i) = sin(in(i));
    });
    Kokkos::fence();

    m_time_stamp += dt;
    m_output.get_header().get_tracking().update_time_stamp(m_time_stamp);
  }

  // Clean up
  void finalize ( ) {}

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const {
    using pack_type = pack::Pack<Real,PackSize>;
    field_repo.template register_field<pack_type>(*m_input_fids.begin());
    field_repo.template register_field<pack_type>(*m_output_fids.begin());
  }

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_input_fids; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_output_fids; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real, device_type>& f) {
    m_input = f;
  }
  void set_computed_field_impl (const Field<      Real, device_type>& f) {
    m_output = f;
  }

  util::TimeStamp           m_time_stamp;

  std::set<FieldIdentifier> m_input_fids;
  std::set<FieldIdentifier> m_output_fids;

  Field<const Real,device_type> m_input;
  Field<Real,device_type>       m_output;

  std::shared_ptr<AbstractGrid> m_grid;

  std::string m_name;

  ParameterList m_params;
  int     m_id;

  Comm    m_comm;
};

} // namespace scream