// sensitivity_output.hpp — partner-side output contract (C++).
//
// Structures and a minimal JSON writer for the black-box response and the
// sensitivity package. Matches the schemas in docs/residual_provider_contract.md
// and docs/output_contract.md. Header-only; no external dependencies.
//
// The response your executable writes must contain at least
// `residual_coefficients` (R^(p), ndof x m, row-major). Everything else is
// optional and stays private unless you choose to share it.

#ifndef RESASM_PARTNER_SENSITIVITY_OUTPUT_HPP
#define RESASM_PARTNER_SENSITIVITY_OUTPUT_HPP

#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace resasm {

struct SensitivityResponse {
  std::string schema = "resasm-partner-response/1";
  std::string status = "ok";                 // "ok" | "error"
  std::size_t ndof = 0;
  std::size_t m = 0;                          // number of seeded parameters
  int order = 1;
  std::vector<double> residual_real;          // ndof   (optional)
  std::vector<double> residual_coefficients;  // ndof*m row-major  (R^(p)) REQUIRED
  std::vector<double> tangent;                // ndof*ndof row-major (optional)
  std::string diagnostics_json = "{}";        // free-form JSON object
  std::string message;
};

namespace detail {
inline void write_array(std::ostream& os, const std::vector<double>& a,
                        std::size_t rows, std::size_t cols) {
  os << "[";
  for (std::size_t i = 0; i < rows; ++i) {
    os << (i ? ",[" : "[");
    for (std::size_t j = 0; j < cols; ++j) {
      if (j) os << ",";
      os << a[i * cols + j];
    }
    os << "]";
  }
  os << "]";
}
inline void write_vec(std::ostream& os, const std::vector<double>& a) {
  os << "[";
  for (std::size_t i = 0; i < a.size(); ++i) { if (i) os << ","; os << a[i]; }
  os << "]";
}
}  // namespace detail

// Write a response.json the kit's BlackBoxRunner can read.
inline bool write_response_json(const std::string& path,
                                const SensitivityResponse& r) {
  std::ofstream os(path);
  if (!os) return false;
  os << "{\n";
  os << "  \"schema\": \"" << r.schema << "\",\n";
  os << "  \"status\": \"" << r.status << "\",\n";
  os << "  \"order\": " << r.order << ",\n";
  os << "  \"residual_coefficients\": ";
  detail::write_array(os, r.residual_coefficients, r.ndof, r.m);
  os << ",\n";
  if (!r.residual_real.empty()) {
    os << "  \"residual_real\": ";
    detail::write_vec(os, r.residual_real);
    os << ",\n";
  }
  if (!r.tangent.empty()) {
    os << "  \"tangent\": ";
    detail::write_array(os, r.tangent, r.ndof, r.ndof);
    os << ",\n";
  }
  os << "  \"diagnostics\": " << r.diagnostics_json << ",\n";
  os << "  \"message\": \"" << r.message << "\"\n";
  os << "}\n";
  return true;
}

// Minimal reader for the request.json the kit writes (fields you need).
// For anything beyond this, use a real JSON library on your side.
struct SensitivityRequest {
  int order = 1;
  std::vector<double> solution;                 // U
  // parameters / seed_directions are name->value maps; parse with your own JSON
  // reader. This header intentionally stays dependency-free.
};

}  // namespace resasm
#endif  // RESASM_PARTNER_SENSITIVITY_OUTPUT_HPP
