// element_residual.hpp — starter templated element residual (edit me).
//
// Keep this generic in `Scalar` so the SAME code runs with double / Dual1 / OTI.
// Replace the body with your (private) element physics. The kit / your driver
// seeds a parameter and reads the imaginary part to obtain a column of R^(1).

#ifndef PARTNER_TEMPLATE_ELEMENT_RESIDUAL_HPP
#define PARTNER_TEMPLATE_ELEMENT_RESIDUAL_HPP

#include <cstddef>

// A cubic axial bar element: N = k * (u1 - u0)^3, residual = [-N, +N].
template <class Scalar>
inline void element_residual(const Scalar* u, const Scalar* params,
                             Scalar* r_e, Scalar* tangent /*may be null*/) {
  const Scalar k = params[0];
  const Scalar d = u[1] - u[0];
  const Scalar N = k * d * d * d;
  r_e[0] = Scalar(0) - N;
  r_e[1] = N;
  if (tangent) {
    // dr/du (real use only): 3 k d^2 * [[1,-1],[-1,1]]
    const Scalar kd = Scalar(3) * k * d * d;
    tangent[0] = kd;  tangent[1] = Scalar(0) - kd;
    tangent[2] = Scalar(0) - kd;  tangent[3] = kd;
  }
}

#endif
