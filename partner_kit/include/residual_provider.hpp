// residual_provider.hpp — partner-side residual provider contract (C++).
//
// Header-only contract a partner can implement to expose a residual evaluator to
// the kit WITHOUT sharing source: they compile their own translation unit and,
// for the black-box path, only ship an executable. The residual is templated on
// the scalar type so the SAME code evaluates with double (real), a first-order
// dual (Dual1), or a future OTI/HYPAD scalar — no rewrite for higher orders.
//
//   template<class Scalar> void eval_element_residual(...);
//   template<class Scalar> void eval_global_residual(...);
//
// Nothing here forces you to reveal the mesh or the material model.

#ifndef RESASM_PARTNER_RESIDUAL_PROVIDER_HPP
#define RESASM_PARTNER_RESIDUAL_PROVIDER_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <cmath>

namespace resasm {

// --------------------------------------------------------------------------
// Minimal first-order dual scalar (reference; swap for your OTI/HYPAD type).
// Any scalar type works if it supports + - * / and integer pow used by your
// residual. Keep your residual written against `Scalar`, not `double`.
// --------------------------------------------------------------------------
struct Dual1 {
  double re{0.0};
  double im{0.0};
  Dual1() = default;
  Dual1(double r) : re(r), im(0.0) {}          // implicit from real
  Dual1(double r, double i) : re(r), im(i) {}
};
inline Dual1 operator+(Dual1 a, Dual1 b){ return {a.re+b.re, a.im+b.im}; }
inline Dual1 operator-(Dual1 a, Dual1 b){ return {a.re-b.re, a.im-b.im}; }
inline Dual1 operator-(Dual1 a){ return {-a.re, -a.im}; }
inline Dual1 operator*(Dual1 a, Dual1 b){ return {a.re*b.re, a.re*b.im + a.im*b.re}; }
inline Dual1 operator/(Dual1 a, Dual1 b){ double c=b.re; return {a.re/c,(a.im*c-a.re*b.im)/(c*c)}; }
inline double real_part(const Dual1& x){ return x.re; }
inline double imag_part(const Dual1& x){ return x.im; }
inline double real_part(double x){ return x; }
inline double imag_part(double){ return 0.0; }
inline Dual1 seed(double v){ return {v, 1.0}; }  // unit first-order perturbation

// --------------------------------------------------------------------------
// Element residual provider. Implement eval_element_residual as a template so it
// can be instantiated for double and Dual1 (and later OTI/HYPAD).
// --------------------------------------------------------------------------
struct ElementDescriptor {
  int element_id{0};
  std::vector<int> dof_map;      // global DOF indices this element touches
};

template <class Scalar>
struct ElementResidualProvider {
  virtual ~ElementResidualProvider() = default;

  // metadata (no proprietary details)
  virtual std::string name() const = 0;
  virtual std::size_t ndof() const = 0;
  virtual std::vector<std::string> parameters() const = 0;
  virtual std::vector<ElementDescriptor> elements() const = 0;

  // r_e has size dof_map.size(); tangent (row-major, optional) is n x n.
  // `params` are indexed to match parameters(); values may carry a Scalar seed.
  virtual void eval_element_residual(
      int element_id,
      const Scalar* element_dofs, std::size_t ndofs,
      const Scalar* params, std::size_t nparams,
      double time, double dtime,
      Scalar* r_e,                         // out: element residual
      Scalar* tangent_optional /*n*n*/,    // out or nullptr
      void*   state_optional) const = 0;   // opaque partner state or nullptr
};

// --------------------------------------------------------------------------
// Global residual provider (whole-model residual vector).
// --------------------------------------------------------------------------
template <class Scalar>
struct GlobalResidualProvider {
  virtual ~GlobalResidualProvider() = default;
  virtual std::string name() const = 0;
  virtual std::size_t ndof() const = 0;
  virtual std::vector<std::string> parameters() const = 0;

  virtual void eval_global_residual(
      const Scalar* U, std::size_t ndof,
      const Scalar* params, std::size_t nparams,
      double time, double dtime,
      Scalar* R /*ndof*/) const = 0;

  // optional tangent; return false if not provided
  virtual bool get_tangent(const double* U, std::size_t ndof,
                           const double* params, std::size_t nparams,
                           double time, double dtime,
                           double* T_rowmajor /*ndof*ndof*/) const { return false; }
  virtual bool apply_tangent(const double* U, std::size_t ndof,
                             const double* params, std::size_t nparams,
                             const double* x, double time, double dtime,
                             double* Tx /*ndof*/) const { return false; }
};

// --------------------------------------------------------------------------
// Convenience: build R^(1) column i by seeding parameter i and reading imag.
// (Order 1; higher orders require an OTI/HYPAD Scalar.)
// --------------------------------------------------------------------------
template <class Provider>
inline void generate_rhs_order1_global(
    const Provider& provider,
    const std::vector<double>& U,
    const std::vector<double>& param_values,
    double time, double dtime,
    std::vector<double>& R1 /*ndof*m, row-major*/) {
  const std::size_t ndof = provider.ndof();
  const std::size_t m = param_values.size();
  R1.assign(ndof * m, 0.0);
  std::vector<Dual1> Ud(U.begin(), U.end());
  std::vector<Dual1> Rd(ndof);
  for (std::size_t j = 0; j < m; ++j) {
    std::vector<Dual1> p(param_values.begin(), param_values.end());
    p[j] = seed(param_values[j]);                    // seed parameter j
    provider.eval_global_residual(Ud.data(), ndof, p.data(), m, time, dtime, Rd.data());
    for (std::size_t i = 0; i < ndof; ++i) R1[i * m + j] = imag_part(Rd[i]);
  }
}

}  // namespace resasm
#endif  // RESASM_PARTNER_RESIDUAL_PROVIDER_HPP
