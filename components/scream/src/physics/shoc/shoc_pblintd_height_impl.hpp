#ifndef SHOC_PBLINTD_HEIGHT_IMPL_HPP
#define SHOC_PBLINTD_HEIGHT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd_height. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd_height(const MemberType&            team,
                                    const Int& nlev,
                                    const Int& npbl,
                                    const uview_1d<const Spack>& z,
                                    const uview_1d<const Spack>& u,
                                    const uview_1d<const Spack>& v,
                                    const Scalar& ustar,
                                    const uview_1d<const Spack>& thv,
                                    const Scalar& thv_ref,
                                    Scalar& pblh,
                                    const uview_1d<Spack>& rino,
                                    bool& check)
{
  const auto s_u = ekat::scalarize(u);
  const auto s_v = ekat::scalarize(v);
  const auto s_rino = ekat::scalarize(rino);
  const auto s_z = ekat::scalarize(z);
  const auto s_thv = ekat::scalarize(thv);

  for (Int k_inv=nlev-2; k_inv >nlev-npbl; --k_inv) {

    if (check==true) {

      const Scalar vvk = ekat::impl::max((s_u(k_inv)-s_u(nlev-1))*(s_u(k_inv)-s_u(nlev-1))
                       + (s_v(k_inv)-s_v(nlev-1))*(s_v(k_inv)-s_v(nlev-1))
                       + (SC::fac*ustar)*(SC::fac*ustar),Scalar(1e-36));

      s_rino(k_inv) = C::gravit*(s_thv(k_inv)-thv_ref)*(s_z(k_inv)-s_z(nlev-1))/(s_thv(nlev-1)*vvk);

      if (s_rino(k_inv) >= SC::ricr) {
        pblh = s_z(k_inv+1) + (SC::ricr - s_rino(k_inv+1))/(s_rino(k_inv)-s_rino(k_inv+1))*
            (s_z(k_inv) - s_z(k_inv+1));
        check = false;
      }

    }


  }
}

} // namespace shoc
} // namespace scream

#endif
