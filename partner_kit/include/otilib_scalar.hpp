// otilib_scalar.hpp — partner-side OTI scalar seam (C++).
//
// Compile your residual ONCE, generic in the scalar type:
//
//     template <class Scalar>
//     void eval_residual(const Scalar* U, const Scalar* params, Scalar* R, ...);
//
// Instantiate it with `double` for the real solve and with an OTILib scalar for
// arbitrary-order sensitivity — no rewrite, no code shared with us.
//
// OTILib (genuine `pyoti` OTI numbers) is an EXTERNAL, GPLv3 dependency you
// install and build yourself from https://github.com/mauriaristi/otilib.git
// (branch master; Windows = WSL only). See docs/otilib_partner_quickstart.md.
// This header does NOT vendor it; it only defines the seam and the extraction
// helpers your driver uses once OTILib's C/C++ headers are available. Define
// RESASM_HAVE_OTILIB and include OTILib before this header to enable the OTI
// typedefs.
//
// Direction convention: exponent multi-index kappa (length m). The true partial
// derivative is (prod_i kappa_i!) * coeff(R, kappa).

#ifndef RESASM_PARTNER_OTILIB_SCALAR_HPP
#define RESASM_PARTNER_OTILIB_SCALAR_HPP

#include <cstddef>
#include <vector>

namespace resasm {

// Real scalar alias (always available).
using Real = double;

#ifdef RESASM_HAVE_OTILIB
// Map these to your OTILib build's real types/functions. The names below are
// placeholders — set them to the genuine OTILib API in ONE place.
//
//   using OtiScalar = oti::number;                 // OTILib scalar type
//   inline OtiScalar oti_seed(double v, int i);    // v + e_i
//   inline double    oti_real(const OtiScalar&);    // real part
//   inline double    oti_coeff(const OtiScalar&, const std::vector<int>& idx);
//
// where `idx` is the basis-index list with multiplicity (kappa -> [1,1,3] etc.).
//
// Provide these in a small otilib_bind.hpp you keep private; then:
//   template <class Scalar> void eval_residual(...);   // your code, unchanged
//   eval_residual<OtiScalar>(...);                     // OTI sensitivity mode
#endif  // RESASM_HAVE_OTILIB

// Helper: exponent multi-index -> basis-index list (kappa=(2,0,1) -> {1,1,3}).
inline std::vector<int> exponents_to_index_list(const std::vector<int>& kappa) {
  std::vector<int> idx;
  for (std::size_t i = 0; i < kappa.size(); ++i)
    for (int r = 0; r < kappa[i]; ++r) idx.push_back(static_cast<int>(i) + 1);
  return idx;
}

// Helper: recovery factor prod_i kappa_i! for a direction.
inline long recovery_factor(const std::vector<int>& kappa) {
  long f = 1;
  for (int k : kappa) { long kf = 1; for (int j = 2; j <= k; ++j) kf *= j; f *= kf; }
  return f;
}

}  // namespace resasm
#endif  // RESASM_PARTNER_OTILIB_SCALAR_HPP
