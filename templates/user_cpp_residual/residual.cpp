// residual.cpp — residual template (Path B: local compiled code).
//
// One binary, two jobs:
//
//   ./residual
//       Self-test. Prints R, the tangent and dR/dk at a known point and
//       cross-checks the analytic parameter derivatives against finite
//       differences of eval_residual. Run it after every edit.
//
//   ./residual --request <req.json> --response <resp.npz>
//       Black-box responder. This is what `resasm run` invokes (see
//       resasm.yml -> residual.command). It reads the framework's request,
//       computes the tangent dR/du and the order-1 right-hand side dR/da_i,
//       and writes the response file.
//
// Write your residual ONCE, generic in the scalar type, and instantiate it with
//   * double            -> real evaluation (what this file does)
//   * an OTILib scalar  -> arbitrary-order sensitivities
// (see ../../partner_kit/include/otilib_scalar.hpp for the OTI seam).
//
// RESPONSE FILE NOTE: the framework hands us a --response path ending in `.npz`.
// A compiled program has no business writing numpy archives, so we strip the
// trailing `.npz` and write `<that>.json` instead. The framework's reader looks
// for exactly that fallback (see resasm_user/providers.py::_read_response).

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

// ===========================================================================
// 1. YOUR MODEL — this is the only section you need to edit.
//
//    R = k u^3 - f        (1-DOF cubic spring)
//    params[0] = k, params[1] = f   (same order as `parameters:` in resasm.yml)
// ===========================================================================

constexpr std::size_t MODEL_NDOF   = 1;   // unknowns          (problem.unknowns)
constexpr std::size_t MODEL_NPARAM = 2;   // parameters, in resasm.yml key order

// The residual itself.
template <class Scalar>
void eval_residual(const Scalar* u, const Scalar* params,
                   std::size_t /*ndof*/, Scalar* R) {
  Scalar k = params[0];
  Scalar f = params[1];
  R[0] = k * (u[0] * u[0] * u[0]) - f;
}

// The tangent dR/du (row-major, ndof x ndof). Used for the sensitivity solve.
template <class Scalar>
void eval_tangent(const Scalar* u, const Scalar* params,
                  std::size_t /*ndof*/, Scalar* T) {
  Scalar k = params[0];
  T[0] = Scalar(3) * k * u[0] * u[0];
}

// dR/d(parameter ip), analytic. `ip` is 0-based, in resasm.yml parameter order.
//   dR/dk = u^3        dR/df = -1
// EDIT THIS WHENEVER YOU EDIT eval_residual. The self-test below finite-
// differences eval_residual and will fail loudly if the two ever drift apart.
template <class Scalar>
void eval_dR_dparam(const Scalar* u, const Scalar* /*params*/,
                    std::size_t /*ndof*/, std::size_t ip, Scalar* dR) {
  if (ip == 0)  dR[0] = u[0] * u[0] * u[0];   // dR/dk
  else          dR[0] = Scalar(-1);           // dR/df
}

// ===========================================================================
// 2. Minimal JSON reading.
//
// The request has a fixed, simple shape, so a tolerant scanner is plenty: we
// only need "u":[...], "parameters":{...}, "order":N and
// "direction_map":{"1":[[...],[...]]}. No escapes, no nesting surprises.
// ===========================================================================

static std::string slurp(const char* path) {
  std::FILE* fh = std::fopen(path, "rb");
  if (!fh) {
    std::fprintf(stderr, "residual: cannot open request file: %s\n", path);
    std::exit(2);
  }
  std::string s;
  char buf[4096];
  std::size_t n;
  while ((n = std::fread(buf, 1, sizeof buf, fh)) > 0) s.append(buf, n);
  std::fclose(fh);
  return s;
}

// Index of the first character of the value of "key", searching from `from`.
static std::size_t value_pos(const std::string& s, const std::string& key,
                             std::size_t from = 0) {
  const std::string tok = "\"" + key + "\"";
  std::size_t p = s.find(tok, from);
  if (p == std::string::npos) return std::string::npos;
  p = s.find(':', p + tok.size());
  if (p == std::string::npos) return std::string::npos;
  return p + 1;
}

// Whitespace and the separators we never care about.
static void skip_sep(const std::string& s, std::size_t& p) {
  while (p < s.size() && (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' ||
                          s[p] == '\r' || s[p] == ','))
    ++p;
}

// One number at/after p; p advances past it.
static double read_number(const std::string& s, std::size_t& p) {
  skip_sep(s, p);
  const char* start = s.c_str() + p;
  char* end = nullptr;
  const double v = std::strtod(start, &end);
  p += static_cast<std::size_t>(end - start);
  return v;
}

// "[a, b, c]" at/after p -> the numbers; p ends just past the closing ']'.
static std::vector<double> read_number_array(const std::string& s, std::size_t& p) {
  std::vector<double> out;
  while (p < s.size() && s[p] != '[') ++p;
  ++p;                                              // past '['
  for (;;) {
    skip_sep(s, p);
    if (p >= s.size() || s[p] == ']') { ++p; break; }
    out.push_back(read_number(s, p));
  }
  return out;
}

// {"k": 2.0, "f": 16.0} -> names + values IN FILE ORDER, which is the
// resasm.yml parameter order (the request's "seed_directions" says the same).
static void read_number_object(const std::string& s, std::size_t& p,
                               std::vector<std::string>& names,
                               std::vector<double>& values) {
  while (p < s.size() && s[p] != '{') ++p;
  ++p;                                              // past '{'
  for (;;) {
    skip_sep(s, p);
    if (p >= s.size() || s[p] == '}') { ++p; break; }
    if (s[p] != '"') { ++p; continue; }
    const std::size_t q = s.find('"', p + 1);
    names.push_back(s.substr(p + 1, q - p - 1));
    p = s.find(':', q) + 1;
    values.push_back(read_number(s, p));
  }
}

// "direction_map": {"1": [[1,0],[0,1]]} -> the exponent vectors for `order`,
// one per column of the order-`order` right-hand side.
static std::vector<std::vector<int>> read_direction_map(const std::string& s,
                                                        int order) {
  std::vector<std::vector<int>> cols;
  const std::size_t dm = value_pos(s, "direction_map");
  if (dm == std::string::npos) return cols;
  std::size_t p = value_pos(s, std::to_string(order), dm);   // the "<order>" key
  if (p == std::string::npos) return cols;
  while (p < s.size() && s[p] != '[') ++p;
  ++p;                                              // past the OUTER '['
  for (;;) {
    skip_sep(s, p);
    if (p >= s.size() || s[p] == ']') { ++p; break; }
    const std::vector<double> v = read_number_array(s, p);   // one exponent vector
    cols.push_back(std::vector<int>(v.begin(), v.end()));
  }
  return cols;
}

// ===========================================================================
// 3. Minimal JSON writing + the black-box responder.
// ===========================================================================

// The framework gives us "<dir>/response.npz"; we write "<dir>/response.json".
static std::string response_json_path(const std::string& resp) {
  std::string out = resp;
  if (out.size() > 4 && out.compare(out.size() - 4, 4, ".npz") == 0)
    out.erase(out.size() - 4);
  if (!(out.size() > 5 && out.compare(out.size() - 5, 5, ".json") == 0))
    out += ".json";
  return out;
}

static void write_matrix(std::FILE* fh, const std::vector<std::vector<double>>& m) {
  std::fputc('[', fh);
  for (std::size_t i = 0; i < m.size(); ++i) {
    if (i) std::fputc(',', fh);
    std::fputc('[', fh);
    for (std::size_t j = 0; j < m[i].size(); ++j) {
      if (j) std::fputc(',', fh);
      std::fprintf(fh, "%.17g", m[i][j]);
    }
    std::fputc(']', fh);
  }
  std::fputc(']', fh);
}

static int respond(const char* req_path, const char* resp_path) {
  const std::string s = slurp(req_path);

  // --- read the request ---------------------------------------------------
  std::size_t p = value_pos(s, "u");
  const std::vector<double> u = read_number_array(s, p);

  std::vector<std::string> pnames;
  std::vector<double> pvals;
  p = value_pos(s, "parameters");
  read_number_object(s, p, pnames, pvals);

  p = value_pos(s, "order");
  const int order = static_cast<int>(read_number(s, p));

  const std::vector<std::vector<int>> dirs = read_direction_map(s, order);

  // --- sanity: this template's model functions are written for one specific
  //     shape. Say so rather than silently returning nonsense. --------------
  if (u.size() != MODEL_NDOF || pvals.size() != MODEL_NPARAM) {
    std::fprintf(stderr,
                 "residual: request has ndof=%zu, nparam=%zu but this model is "
                 "written for ndof=%zu, nparam=%zu. Update MODEL_NDOF / "
                 "MODEL_NPARAM and the model functions.\n",
                 u.size(), pvals.size(), MODEL_NDOF, MODEL_NPARAM);
    return 2;
  }
  if (order != 1) {
    std::fprintf(stderr,
                 "residual: this template implements order 1 only (got order=%d).\n"
                 "For order >= 2 you must also consume \"u_star_coefficients\" to "
                 "rebuild u* — or use the OTILib scalar type.\n", order);
    return 2;
  }
  if (dirs.empty()) {
    std::fprintf(stderr, "residual: request has no direction_map[\"%d\"].\n", order);
    return 2;
  }

  // --- compute ------------------------------------------------------------
  // tangent dR/du (ndof x ndof)
  std::vector<double> Tflat(MODEL_NDOF * MODEL_NDOF, 0.0);
  eval_tangent<double>(u.data(), pvals.data(), MODEL_NDOF, Tflat.data());
  std::vector<std::vector<double>> T(MODEL_NDOF, std::vector<double>(MODEL_NDOF));
  for (std::size_t i = 0; i < MODEL_NDOF; ++i)
    for (std::size_t j = 0; j < MODEL_NDOF; ++j)
      T[i][j] = Tflat[i * MODEL_NDOF + j];

  // order-1 RHS: one column per direction; R1[:,c] = dR/da_i, where i is the
  // position of the 1 in that direction's exponent vector.
  std::vector<std::vector<double>> R1(MODEL_NDOF,
                                      std::vector<double>(dirs.size(), 0.0));
  std::vector<double> col(MODEL_NDOF, 0.0);
  for (std::size_t c = 0; c < dirs.size(); ++c) {
    std::size_t ip = 0;
    for (std::size_t e = 0; e < dirs[c].size(); ++e)
      if (dirs[c][e] != 0) { ip = e; break; }
    eval_dR_dparam<double>(u.data(), pvals.data(), MODEL_NDOF, ip, col.data());
    for (std::size_t i = 0; i < MODEL_NDOF; ++i) R1[i][c] = col[i];
  }

  // --- write the response -------------------------------------------------
  const std::string out = response_json_path(resp_path);
  std::FILE* fh = std::fopen(out.c_str(), "w");
  if (!fh) {
    std::fprintf(stderr, "residual: cannot write response file: %s\n", out.c_str());
    return 2;
  }
  std::fprintf(fh, "{\"residual_coefficients_by_order\": {\"%d\": ", order);
  write_matrix(fh, R1);
  std::fprintf(fh, "}, \"tangent\": ");
  write_matrix(fh, T);
  std::fprintf(fh,
               ", \"diagnostics\": {\"solver\": \"template_cpp\", "
               "\"method\": \"analytic\", \"ndof\": %zu, \"order\": %d}}\n",
               MODEL_NDOF, order);
  std::fclose(fh);
  return 0;
}

// ===========================================================================
// 4. Self-test (no arguments).
// ===========================================================================

static int self_test() {
  double u[MODEL_NDOF]   = {2.0};
  double p[MODEL_NPARAM] = {2.0, 16.0};        // k, f
  double R[MODEL_NDOF], T[MODEL_NDOF * MODEL_NDOF];

  eval_residual<double>(u, p, MODEL_NDOF, R);
  eval_tangent<double>(u, p, MODEL_NDOF, T);
  std::printf("R = %.6f  (expect 0)\n", R[0]);
  std::printf("T = dR/du = %.6f  (expect 24)\n", T[0]);

  // dR/dk = u^3 = 8 ; du/dk = -(dR/dk)/T = -8/24 = -0.3333
  double dk[MODEL_NDOF];
  eval_dR_dparam<double>(u, p, MODEL_NDOF, 0, dk);
  std::printf("dR/dk = %.6f  du/dk = %.6f\n", dk[0], -dk[0] / T[0]);

  // Guard against the analytic dR/da drifting away from eval_residual: central
  // finite differences must reproduce eval_dR_dparam.
  int bad = 0;
  for (std::size_t ip = 0; ip < MODEL_NPARAM; ++ip) {
    const double h = 1e-6 * std::fmax(1.0, std::fabs(p[ip]));
    double pp[MODEL_NPARAM], pm[MODEL_NPARAM], Rp[MODEL_NDOF], Rm[MODEL_NDOF],
           an[MODEL_NDOF];
    std::memcpy(pp, p, sizeof p);
    std::memcpy(pm, p, sizeof p);
    pp[ip] += h;
    pm[ip] -= h;
    eval_residual<double>(u, pp, MODEL_NDOF, Rp);
    eval_residual<double>(u, pm, MODEL_NDOF, Rm);
    eval_dR_dparam<double>(u, p, MODEL_NDOF, ip, an);
    for (std::size_t i = 0; i < MODEL_NDOF; ++i) {
      const double fd = (Rp[i] - Rm[i]) / (2.0 * h);
      if (std::fabs(fd - an[i]) > 1e-6 * std::fmax(1.0, std::fabs(fd))) {
        std::printf("MISMATCH dR/da[%zu][%zu]: analytic %.9g vs FD %.9g\n",
                    ip, i, an[i], fd);
        bad = 1;
      }
    }
  }
  std::printf("analytic dR/da vs finite differences: %s\n", bad ? "FAIL" : "ok");
  return bad;
}

// ===========================================================================
int main(int argc, char** argv) {
  const char* req = nullptr;
  const char* resp = nullptr;
  for (int i = 1; i < argc - 1; ++i) {
    if (std::strcmp(argv[i], "--request") == 0)  req  = argv[++i];
    else if (std::strcmp(argv[i], "--response") == 0) resp = argv[++i];
  }
  if (req && resp) return respond(req, resp);        // black-box responder
  if (req || resp) {
    std::fprintf(stderr, "usage: residual [--request <in.json> --response <out.npz>]\n");
    return 2;
  }
  return self_test();                                // no args -> self-test
}
