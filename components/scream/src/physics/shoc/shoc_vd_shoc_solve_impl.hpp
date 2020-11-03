#ifndef SHOC_VD_SHOC_SOLVE_HPP
#define SHOC_VD_SHOC_SOLVE_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template <typename TridiagArray>
Kokkos::View<typename TridiagArray::value_type**, Kokkos::LayoutRight>
KOKKOS_INLINE_FUNCTION
get_diags (const TridiagArray& A, const int diag_idx) {
  return Kokkos::View<typename TridiagArray::value_type**, Kokkos::LayoutRight>(
    &A.impl_map().reference(diag_idx, 0, 0),
    A.extent_int(1), A.extent_int(2));
}

template <typename Scalar>
using TridiagArray = Kokkos::View<Scalar***, Kokkos::LayoutRight>;
template <typename Scalar>
using DataArray = Kokkos::View<Scalar**, Kokkos::LayoutRight>;



template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::vd_shoc_solve(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<Spack>& ca,
  const uview_1d<Spack>& cc,
  const uview_1d<const Spack>& denom,
  const uview_1d<const Spack>& ze,
  const uview_1d<Spack>&       diag,
  const uview_1d<Spack>&       var)
{
  using ekat::scalarize;
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;
  using Kokkos::subview;
  using Kokkos::ALL;

  TridiagArray<Scalar> A("A", 3, nlev, 1);
  DataArray<Scalar> B("B", nlev, 1);

  // Fill data
  const auto Am = create_mirror_view(A);
  const auto Bm = create_mirror_view(B);
  {
    //const auto As = scalarize(Am);
    const auto dl = subview(A, 0, ALL(), ALL());
    const auto d  = subview(A, 1, ALL(), ALL());
    const auto du = subview(A, 2, ALL(), ALL());

    for (int i=0; i<1; ++i) {
      for (int k=0; k<nlev; ++k) {
        dl(k,i) = scalarize(cc)(k);
        du(k,i) = scalarize(ca)(k);
        d(k,i) = 100;
      }
    }

    //const auto Bs = scalarize(Bm);
    for (int i=0; i<1; ++i) {
      for (int k=0; k<nlev; ++k) {
        B(k,i) = scalarize(var)(k);
      }
    }
  }

  {
    const auto dl = get_diags(A, 0);
    const auto d  = get_diags(A, 1);
    const auto du = get_diags(A, 2);

    std::cout << std::endl << "DL:  ";
    for (int k=0; k<nlev; ++k) std::cout << dl(k,0) << "   ";
    std::cout << std::endl;
    std::cout << std::endl << "DU:  ";
    for (int k=0; k<nlev; ++k) std::cout << du(k,0) << "   ";
    std::cout << std::endl;
    std::cout << std::endl << "D:   ";
    for (int k=0; k<nlev; ++k) std::cout << d(k,0) << "   ";
    std::cout << std::endl;
    std::cout << std::endl << "RHS: ";
    for (int k=0; k<nlev; ++k) std::cout << B(k,0) << "   ";
    std::cout << std::endl;

    ekat::bfb(team, dl, d, du, B);

    std::cout << std::endl << std::endl << "SOLUTION: ";
    for (int k=0; k<nlev; ++k) std::cout << B(k,0) << "   ";
    std::cout << std::endl;
  }

  for (int i=0; i<1; ++i) {
    for (int k=0; k<nlev; ++k) {
      const int p = k/Spack::n;
      const int p_indx = k%Spack::n;

      var(p)[p_indx] = B(k,0);
    }
  }
}

} // namespace shoc
} // namespace scream

#endif
