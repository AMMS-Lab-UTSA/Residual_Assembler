// residual.hpp — templated residual for real + OTILib sensitivity (edit me).
//
// Keep it generic in Scalar. `params` are seeded (a_i + e_i) in OTI mode; the
// coefficients of R along each order-p direction are R^(p). See the README and
// ../../include/otilib_scalar.hpp.

#ifndef PARTNER_TEMPLATE_OTILIB_RESIDUAL_HPP
#define PARTNER_TEMPLATE_OTILIB_RESIDUAL_HPP

#include <cstddef>

// Example: 1-DOF cubic spring R = k u^3 - f.  params[0]=k, params[1]=f.
template <class Scalar>
inline void eval_residual(const Scalar* U, const Scalar* params,
                          std::size_t /*ndof*/, Scalar* R) {
  Scalar u = U[0];
  Scalar k = params[0];
  Scalar f = params[1];
  R[0] = k * (u * u * u) - f;   // ordinary arithmetic -> double or OTILib scalar
}

// Optional real tangent (dR/dU) for the sensitivity solve; OTILib can also
// provide this by seeding the DOFs instead of the parameters.
template <class Scalar>
inline void eval_tangent(const Scalar* U, const Scalar* params,
                         std::size_t /*ndof*/, Scalar* T /*ndof*ndof*/) {
  Scalar u = U[0];
  Scalar k = params[0];
  T[0] = Scalar(3) * k * u * u;   // 3 k u^2
}

#endif
