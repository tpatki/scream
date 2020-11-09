#ifndef SHOC_VD_SHOC_DECOMP_HPP
#define SHOC_VD_SHOC_DECOMP_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::vd_shoc_decomp(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const uview_1d<const Spack>& kv_term,
  const uview_1d<const Spack>& tmpi,
  const uview_1d<const Spack>& rdp_zt,
  const Scalar&                dtime,
  const Scalar&                flux,
  const uview_1d<Spack>&       ca,
  const uview_1d<Spack>&       cc,
  const uview_1d<Spack>&       denom,
  const uview_1d<Spack>&       ze)
{

}

} // namespace shoc
} // namespace scream

#endif
