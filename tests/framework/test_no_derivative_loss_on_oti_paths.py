"""No operation that loses derivatives enters the OTI chain unreviewed.

A hypercomplex residual keeps its derivatives only while every operation on
it is one OTILib overloads. The ways to lose them are few and recognisable in
source: ``float()`` and ``.real`` (both silently return the real part of a
Dual1 number; ``.real`` also of an OTILib number), ``dtype=float`` and
``astype(float)``, arrays allocated real (``np.zeros``/``np.eye`` without
``dtype=object``), ``math`` functions, ``numpy.linalg``/SciPy, object-array
conversion and serialisation (``json``, ``np.save``, ``pickle``) before the
coefficients are extracted.

This scan walks the syntax tree of every module through which a live OTI
number passes, from seeding (``oti_rhs_provider``, ``finite_strain_sensitivity``,
``resasm_user.oti_global``) through the assembler, the OTI-differentiable
element backends, the finite-strain element, its material and the C3D8
kernels, to the real-only consumers named by the requirement
(``field_sensitivity``, ``replay/j2_history``). Every such operation found
there has been reviewed and is listed in ``REVIEWED`` with the reason it
cannot lose a derivative: it acts on real geometry, on the real tangent, on
coefficients already extracted, or on inputs refused unless real. A new one,
anywhere in these modules, fails this test until it is reviewed and listed;
a listed one that disappears fails it too, so the list stays exact.
"""
import ast
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]

#: The modules a live OTI number flows through, and what each does with it.
MODULES = {
    "residual_core/core/oti_rhs_provider.py": "seeds parameters, evaluates R*, extracts R^(p)",
    "residual_core/formulations/finite_strain_sensitivity.py": "seeds finite-strain PROPS",
    "resasm_user/oti_global.py": "the user-residual order loop",
    "residual_core/core/assembler.py": "scatters OTI element residuals",
    "residual_core/formulations/nonlinear_spring1.py": "OTI-differentiable element",
    "residual_core/formulations/nonlinear_bar1.py": "OTI-differentiable element",
    "residual_core/formulations/solid_c3d8_finite_strain.py": "carries OTI stress",
    "residual_core/materials/neo_hookean.py": "evaluated with OTI PROPS",
    "residual_core/formulations/c3d8_kernel.py": "integrates OTI stress",
    "residual_core/formulations/c3d8_sensitivity.py": "assembles dR/dp, OTI-typed or real",
    "residual_core/core/field_sensitivity.py": "real consumer: refuses live OTI",
    "residual_core/replay/j2_history.py": "real consumer: refuses live OTI",
}

GEOMETRY = "real geometry: coordinates, shape functions and Jacobians are never seeded"
EXTRACTED = "acts on extracted real coefficients and the real tangent, after extraction"
SELF_CHECK = "real verification helper of the kernel; never called with OTI values"
REFUSES = "inputs are validated real (live OTI is refused with an error) before this"

#: (module, function) -> ({operation: count}, why it cannot lose a derivative)
REVIEWED = {
    ("residual_core/core/oti_rhs_provider.py", "OtiLibRHSProvider.evaluate_rhs"): (
        {"float()": 1, "real-only allocation": 3, "linear algebra": 1},
        "float() reads the REAL converged U before it is lifted to OTI; the free "
        "mask, R^(p) and U^(p) are real arrays of extracted coefficients, solved "
        "with the real tangent"),
    ("residual_core/core/oti_rhs_provider.py", "OtiLibRHSProvider._real_parameter_values"): (
        {"float()": 1}, "reads the real parameter value that is then seeded"),
    ("resasm_user/oti_global.py", "solve_python"): (
        {"float()": 2, "real-only allocation": 2, "linear algebra": 1},
        "float() reads the real parameter values and real u before seeding; R^(p) "
        "and U^(p) hold extracted coefficients solved with the real tangent"),
    ("resasm_user/oti_global.py", "solve_executable"): (
        {"real-only allocation": 1, "linear algebra": 1},
        "the executable returns real coefficients; solved with the real tangent"),
    ("residual_core/core/assembler.py", "Assembler.assemble"): (
        {"dtype=float": 1, "real-only allocation": 2},
        "U is the real solution (geometry is not seeded); R and K start real and "
        "are promoted to object arrays as soon as an element returns OTI values"),
    ("residual_core/core/assembler.py", "Assembler.split"): (
        {"real-only allocation": 1}, "an empty reaction vector when nothing is prescribed"),
    ("residual_core/formulations/nonlinear_spring1.py", "_real_of"): (
        {"float()": 2, ".real": 1},
        "diagnostics only: the displacement echoed in the element's info dict, "
        "never part of the residual"),
    ("residual_core/formulations/nonlinear_bar1.py", "NonlinearBar1.eval_element"): (
        {"float()": 2},
        "tangent branch only; the OTI residual is evaluated with compute_tangent=False"),
    ("residual_core/formulations/solid_c3d8_finite_strain.py", "_F_at"): (
        {"real-only allocation": 2}, GEOMETRY),
    ("residual_core/formulations/solid_c3d8_finite_strain.py",
     "SolidC3D8FiniteStrain.eval_element"): (
        {"real-only allocation": 2, "linear algebra": 1},
        "the previous displacement, the Jacobian check and the empty history are "
        "real geometry and state; the OTI stress passes through untouched"),
    ("residual_core/materials/neo_hookean.py", "CompressibleNeoHookean.evaluate"): (
        {"float()": 1, ".real": 1, "linear algebra": 2, "real-only allocation": 3},
        "float(.real) only validates the primal of mu and lambda; det(F) and the "
        "identity and unit-strain tensors are real kinematics; stress and tangent "
        "are built from the OTI constants by overloaded arithmetic"),
    ("residual_core/formulations/c3d8_kernel.py", "<module>"): (
        {"dtype=float": 2}, "constant node and coordinate tables"),
    ("residual_core/formulations/c3d8_kernel.py", "_abaqus_c3d8_gauss"): (
        {"dtype=float": 2}, "quadrature points and weights"),
    ("residual_core/formulations/c3d8_kernel.py", "shape_functions"): (
        {"dtype=float": 1}, GEOMETRY),
    ("residual_core/formulations/c3d8_kernel.py", "shape_grad_natural"): (
        {"dtype=float": 2}, GEOMETRY),
    ("residual_core/formulations/c3d8_kernel.py", "_b_from_coords"): (
        {"dtype=float": 2, "linear algebra": 2}, GEOMETRY),
    ("residual_core/formulations/c3d8_kernel.py", "_as_ue_2d"): (
        {"dtype=float": 1}, "real displacements: geometry is not seeded"),
    ("residual_core/formulations/c3d8_kernel.py", "element_internal_force_small_strain"): (
        {"dtype=float": 1}, GEOMETRY + "; the force takes the stress's own dtype"),
    ("residual_core/formulations/c3d8_kernel.py", "element_internal_force_finite_strain"): (
        {"dtype=float": 1}, GEOMETRY + "; the force takes the stress's own dtype"),
    ("residual_core/formulations/c3d8_kernel.py", "kirchhoff_jaumann_to_spatial"): (
        {"real-only allocation": 1},
        "the unit strain direction is a real constant; the corrected modulus takes "
        "the dtype of DDSDDE and stress, OTI when the PROPS are seeded"),
    ("residual_core/formulations/c3d8_kernel.py", "element_tangent"): (
        {"dtype=float": 1, "real-only allocation": 2},
        GEOMETRY + " (with the zero default displacement and the identity of the "
        "kron product); K takes the dtype of the modulus and stress"),
    ("residual_core/formulations/c3d8_kernel.py", "force_tangent_fixed_sigma"): (
        {"dtype=float": 3}, "a real tangent at a frozen real stress"),
    ("residual_core/formulations/c3d8_kernel.py", "assemble_global_internal_force"): (
        {"dtype=float": 5}, "real global force of the stress-driven check path"),
    ("residual_core/formulations/c3d8_kernel.py", "assemble_global_internal_force.get_u"): (
        {"dtype=float": 1}, "real displacements of the stress-driven check path"),
    ("residual_core/formulations/c3d8_kernel.py", "surface_traction_nodal_forces"): (
        {"dtype=float": 2}, SELF_CHECK),
    ("residual_core/formulations/c3d8_kernel.py", "linear_stress_target"): (
        {"dtype=float": 2, "real-only allocation": 1}, SELF_CHECK),
    ("residual_core/formulations/c3d8_kernel.py", "linear_stress_target.sig_tensor"): (
        {"dtype=float": 1}, SELF_CHECK),
    ("residual_core/formulations/c3d8_kernel.py", "isotropic_D"): (
        {"dtype=float": 1}, "real elastic moduli of the self-checks"),
    ("residual_core/formulations/c3d8_kernel.py", "_check_divergence_theorem"): (
        {"linear algebra": 3}, SELF_CHECK),
    ("residual_core/formulations/c3d8_kernel.py", "_tangent_fd"): (
        {"real-only allocation": 3, "linear algebra": 3}, SELF_CHECK),
    ("residual_core/formulations/c3d8_kernel.py", "_check_finite_strain_fd"): (
        {"real-only allocation": 2, "linear algebra": 3}, SELF_CHECK),
    ("residual_core/formulations/c3d8_kernel.py", "_check_linear_stress"): (
        {"linear algebra": 2}, SELF_CHECK),
    ("residual_core/formulations/c3d8_sensitivity.py", "element_dR_dp"): (
        {"dtype=float": 2}, GEOMETRY + "; dR/dp takes the stress derivative's own dtype"),
    ("residual_core/formulations/c3d8_sensitivity.py", "assemble_dR_dp"): (
        {"dtype=float": 6},
        GEOMETRY + "; the global vector is promoted to the element dtype before scatter"),
    ("residual_core/formulations/c3d8_sensitivity.py", "SensitivityResult.as_dict"): (
        {"float()": 3}, "reports norms of the solved real du/dp"),
    ("residual_core/formulations/c3d8_sensitivity.py", "solve_du_dp"): (
        {"dtype=float": 2, "real-only allocation": 2, "linear algebra": 3, "float()": 2},
        EXTRACTED),
    ("residual_core/core/field_sensitivity.py", "solve_field_sensitivities"): (
        {"dtype=float": 3, "linear algebra": 1}, REFUSES),
    ("residual_core/replay/j2_history.py", "real_array"): (
        {"astype(float)": 1}, "the refusal itself: casts only after checking the dtype is real"),
    ("residual_core/replay/j2_history.py", "step_j2"): (
        {"real-only allocation": 3}, "zero increments of the real elastic-predictor probe"),
}

_ALLOCATORS = {"zeros", "ones", "empty", "full", "zeros_like", "ones_like",
               "empty_like", "full_like", "eye", "identity"}
_SERIALISERS = {("json", "dump"), ("json", "dumps"), ("np", "save"), ("np", "savez"),
                ("np", "savez_compressed"), ("np", "savetxt"), ("numpy", "save"),
                ("pickle", "dump"), ("pickle", "dumps")}
_FLOAT_TYPES = {"float", "np.float64", "numpy.float64", "np.double", "np.float32",
                "np.float_"}


def _dotted(node):
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
        return ".".join(reversed(parts))
    return None


def _is_float_type(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value.startswith("float")
    return _dotted(node) in _FLOAT_TYPES


class _Scanner(ast.NodeVisitor):
    """Collect (function qualname, operation, line) for derivative-losing operations."""

    def __init__(self):
        self.scope, self.found = [], []

    def _add(self, node, operation):
        self.found.append((".".join(self.scope) or "<module>", operation, node.lineno))

    def visit_FunctionDef(self, node):
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    visit_AsyncFunctionDef = visit_FunctionDef
    visit_ClassDef = visit_FunctionDef

    def visit_Attribute(self, node):
        if node.attr == "real" and isinstance(node.ctx, ast.Load):
            self._add(node, ".real")
        self.generic_visit(node)

    def visit_Call(self, node):
        name = _dotted(node.func) or ""
        head, _, tail = name.rpartition(".")
        dtype = [keyword.value for keyword in node.keywords if keyword.arg == "dtype"]
        if name in ("float", "complex"):
            self._add(node, "float()")
        if (name == "getattr" and len(node.args) >= 2
                and isinstance(node.args[1], ast.Constant) and node.args[1].value == "real"):
            self._add(node, ".real")
        if name.startswith("math."):
            self._add(node, "math.")
        if name.startswith(("np.linalg.", "numpy.linalg.", "scipy.")):
            self._add(node, "linear algebra")
        if (isinstance(node.func, ast.Attribute) and node.func.attr == "astype"
                and node.args and _is_float_type(node.args[0])):
            self._add(node, "astype(float)")
        if dtype and _is_float_type(dtype[0]):
            self._add(node, "dtype=float")
        if dtype and _dotted(dtype[0]) == "object":
            self._add(node, "object array")
        if head.split(".")[0] in ("np", "numpy") and tail in _ALLOCATORS and not dtype:
            self._add(node, "real-only allocation")
        if tuple(name.split(".")[-2:]) in _SERIALISERS:
            self._add(node, "serialisation")
        self.generic_visit(node)


def scan(source):
    scanner = _Scanner()
    scanner.visit(ast.parse(source))
    return scanner.found


def found_in_modules():
    """{(module, function): Counter(operation)} and the lines behind each."""
    counts, lines = defaultdict(Counter), defaultdict(list)
    for module in MODULES:
        for function, operation, line in scan((ROOT / module).read_text()):
            counts[(module, function)][operation] += 1
            lines[(module, function)].append("%s:%d %s" % (module, line, operation))
    return counts, lines


def test_the_scanner_recognises_every_derivative_losing_operation():
    source = '''
import json, math, pickle
import numpy as np
def chain(x, xs):
    a = float(x)
    b = x.real
    c = getattr(x, "real", x)
    d = np.asarray(xs, dtype=float)
    e = np.array(xs).astype(np.float64)
    f = np.zeros(3)
    g = math.exp(x)
    h = np.linalg.solve(xs, xs)
    i = np.array(xs, dtype=object)
    json.dumps(xs); np.save("f", xs); pickle.dumps(xs)
    # what keeps derivatives is not reported
    kept = np.zeros(3, dtype=np.result_type(xs.dtype, float)) + x * x
    return np.asarray(xs).astype(object) + kept
'''
    operations = Counter(operation for _, operation, _ in scan(source))
    assert operations == {"float()": 1, ".real": 2, "dtype=float": 1, "astype(float)": 1,
                          "real-only allocation": 1, "math.": 1, "linear algebra": 1,
                          "object array": 1, "serialisation": 3}


def test_every_derivative_losing_operation_on_the_oti_chain_is_reviewed():
    counts, lines = found_in_modules()
    unreviewed = []
    for key, found in sorted(counts.items()):
        reviewed = REVIEWED.get(key, ({}, None))[0]
        if dict(found) != reviewed:
            unreviewed.append("%s in %s: found %s, reviewed %s\n    %s" % (
                key[1], key[0], dict(found), reviewed, "\n    ".join(lines[key])))
    assert not unreviewed, (
        "derivative-losing operations on the OTI chain that nobody has reviewed "
        "(prove each cannot see a live OTI number, then list it in REVIEWED):\n"
        + "\n".join(unreviewed))


def test_the_reviewed_list_names_only_what_is_there():
    counts, _ = found_in_modules()
    stale = [key for key in REVIEWED if key not in counts]
    assert not stale, "reviewed entries with nothing left to review: %s" % stale
    assert {module for module, _ in REVIEWED} <= set(MODULES)
    for (module, function), (_, reason) in REVIEWED.items():
        assert reason, (module, function)


def test_the_modules_scanned_are_where_oti_numbers_flow():
    """The chain is not hypothetical: seeding and the element and material
    entry points named above exist and are what the providers call."""
    for module in MODULES:
        assert (ROOT / module).is_file(), module
    source = (ROOT / "residual_core/core/oti_rhs_provider.py").read_text()
    assert "ctx.seed(" in source and "form.eval_element(" in source
    source = (ROOT / "residual_core/formulations/finite_strain_sensitivity.py").read_text()
    assert "context.seed(" in source and "Assembler(" in source


def test_a_genuine_oti_number_refuses_the_loud_operations():
    """Why ``.real`` and ``float()`` are the silent risks: with genuine
    OTILib the array casts, math, linear algebra and real scatter raise
    instead of dropping the derivative, while ``.real`` does not."""
    from residual_core.algebra import otilib_adapter
    if not otilib_adapter.otilib_available():
        pytest.skip("OTILib not installed")
    number = otilib_adapter.OtiContext(1, 1).seed(2.0, 1)
    for operation in (lambda: float(number),
                      lambda: np.asarray([number], dtype=float),
                      lambda: np.linalg.det(np.array([[number]], dtype=object)),
                      lambda: np.add.at(np.zeros(1), [0], np.array([number], dtype=object))):
        with pytest.raises((TypeError, ValueError)):
            operation()
    assert number.real == 2.0
