// main.cpp — C++ black-box residual solver skeleton (edit the residual).
//
// Speaks the kit's contract:  your_solver --input request.json --output response.json
// Reads the requested solution + parameter values + seeded directions, evaluates a
// TEMPLATED global residual with resasm::Dual1, and writes R^(1) to response.json.
//
// NOTE: the JSON reading below is a minimal illustrative extractor (numbers only)
// so this file is dependency-free. Replace with a real JSON library for production.
//
//   c++ -std=c++17 -I../../include main.cpp -o your_solver

#include "residual_provider.hpp"
#include "sensitivity_output.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using resasm::Dual1;

// ---- your private residual: keep it generic in Scalar ----------------------
template <class Scalar>
void eval_global_residual(const Scalar* U, std::size_t /*ndof*/,
                          const Scalar* params, std::size_t /*np*/,
                          Scalar* R) {
  // Example: 1-DOF cubic spring R = k u^3 - f
  Scalar u = U[0];
  Scalar k = params[0];
  Scalar f = params[1];
  R[0] = k * (u * u * u) - f;
}

// ---- tiny helpers to pull numbers out of the request JSON (illustrative) ---
static std::string slurp(const char* path) {
  std::ifstream is(path);
  std::stringstream ss; ss << is.rdbuf(); return ss.str();
}
// read the numeric array following a "key": [ ... ]
static std::vector<double> read_array(const std::string& s, const std::string& key) {
  std::vector<double> out;
  auto k = s.find("\"" + key + "\"");
  if (k == std::string::npos) return out;
  auto lb = s.find('[', k); auto rb = s.find(']', lb);
  if (lb == std::string::npos || rb == std::string::npos) return out;
  std::string body = s.substr(lb + 1, rb - lb - 1);
  std::stringstream ss(body); std::string tok;
  while (std::getline(ss, tok, ',')) { try { out.push_back(std::stod(tok)); } catch (...) {} }
  return out;
}
// read a single number following a "key": value
static double read_scalar(const std::string& s, const std::string& key, double dflt) {
  auto k = s.find("\"" + key + "\"");
  if (k == std::string::npos) return dflt;
  auto c = s.find(':', k);
  try { return std::stod(s.substr(c + 1)); } catch (...) { return dflt; }
}

int main(int argc, char** argv) {
  const char* in = nullptr; const char* out = nullptr;
  for (int i = 1; i + 1 < argc; ++i) {
    if (!std::strcmp(argv[i], "--input")) in = argv[i + 1];
    if (!std::strcmp(argv[i], "--output")) out = argv[i + 1];
  }
  if (!in || !out) { std::fprintf(stderr, "need --input and --output\n"); return 2; }

  std::string req = slurp(in);
  std::vector<double> U = read_array(req, "solution");
  // parameters come as an object; for this skeleton read k and f by name
  double k = read_scalar(req, "k", 1.0);
  double f = read_scalar(req, "f", 0.0);

  const std::size_t ndof = U.size() ? U.size() : 1;

  // --- seed the design parameter k (order 1) and evaluate ---
  std::vector<Dual1> Ud(U.begin(), U.end());
  if (Ud.empty()) Ud.push_back(Dual1(0.0));
  std::vector<Dual1> p = { resasm::seed(k), Dual1(f) };   // seed k along eps
  std::vector<Dual1> Rd(ndof);
  eval_global_residual<Dual1>(Ud.data(), ndof, p.data(), p.size(), Rd.data());

  resasm::SensitivityResponse resp;
  resp.ndof = ndof; resp.m = 1; resp.order = 1;
  resp.residual_real.resize(ndof);
  resp.residual_coefficients.resize(ndof * 1);
  for (std::size_t i = 0; i < ndof; ++i) {
    resp.residual_real[i] = resasm::real_part(Rd[i]);
    resp.residual_coefficients[i] = resasm::imag_part(Rd[i]);   // dR/dk
  }
  resp.diagnostics_json = "{\"method\":\"cpp-dual1\"}";
  resasm::write_response_json(out, resp);
  return 0;
}
