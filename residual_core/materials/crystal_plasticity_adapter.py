"""Crystal-plasticity specialization of the generic UMAT bridge.

Configures :class:`~residual_core.materials.umat_adapter.UmatAdapter` for the
Oxford / Grilli crystal-plasticity UMAT
(``sources/permissive/ngrilli_Oxford_Crystal_Plasticity/umat.for``): 125 state
dependent variables, an 11-constant PROPS block, and the STATEV map inferred
from ``umat.for`` / ``kmat.f``.

PROPS layout (11 constants; 0-based here, 1-based in the .inp)
    [0]    crystal type   (0=HCP, 1/2=BCC, 3=Carbide, 4=Olivine, 5=Orthorhombic)
    [1:10] rotation matrix crystal->sample, order R11,R12,R13,R21,R22,R23,R31,R32,R33
    [10]   grain index

STATEV (SDV) map used here (0-based slices; see umat_adapter_fortran/README.md)
    [0:9]    rotation matrix gmatinv (crystal->sample)      (SDV 1-9)
    [34]     cumulative plastic slip                        (SDV 35)
    [47:53]  Cauchy stress s11,s22,s33,s12,s13,s23          (SDV 48-53)
    [80:89]  plastic deformation gradient Fp (row-major)    (SDV 81-89)

HONESTY / CAPABILITY
    The ACTUAL crystal-plasticity stress requires the compiled REAL UMAT, which
    only builds with Intel ifort + Abaqus (MKL); gfortran cannot compile the
    unmodified sources (see umat_adapter_fortran/README.md 'Build status').
    OFFLINE (no Abaqus) this class can only run through backend='python' with an
    injected MOCK UMAT -- which is NOT crystal plasticity, only a wiring proxy.
    With backend='fortran' and no real binary it raises a clear error and never
    fabricates a crystal-plasticity stress.
"""

from __future__ import annotations

from typing import Any, Callable, Optional

import numpy as np

from .base import MaterialBinding
from .umat_adapter import UmatAdapter


class CrystalPlasticityAdapter(UmatAdapter):
    """Grilli Oxford CP UMAT behind the Material interface."""

    name = "crystal_plasticity"

    # --- fixed model sizes ------------------------------------------------
    N_STATE_VARS = 125
    N_PROPS = 11

    # --- PROPS layout (0-based) ------------------------------------------
    IDX_CRYSTAL_TYPE = 0
    SLICE_ROTATION_PROPS = slice(1, 10)     # PROPS 2-10
    IDX_GRAIN_ID = 10                       # PROPS 11

    # --- STATEV map (0-based) --------------------------------------------
    SLICE_ROTATION = slice(0, 9)            # SDV 1-9
    IDX_CUM_SLIP = 34                       # SDV 35
    SLICE_STRESS = slice(47, 53)            # SDV 48-53
    SLICE_FP = slice(80, 89)                # SDV 81-89
    IDX_SSD = 53                            # SDV 54 (sessile SSD density)

    CRYSTAL_TYPES = {0: "HCP", 1: "BCC", 2: "BCC", 3: "Carbide",
                     4: "Olivine", 5: "Orthorhombic"}

    def __init__(self,
                 backend: str = "fortran",
                 umat_fn: Optional[Callable[..., Any]] = None,
                 driver: Optional[str] = None,
                 allow_mock_driver: bool = False,
                 **kwargs):
        # allow_mock_driver defaults False: the isotropic mock binary is NOT
        # crystal plasticity, so fortran backend must find the real ifort+Abaqus
        # binary or raise.
        super().__init__(
            n_state_vars=self.N_STATE_VARS,
            backend=backend,
            umat_fn=umat_fn,
            driver=driver,
            cmname=kwargs.pop("cmname", "CPMATERIAL"),
            rotation_prop_offset=1,
            rotation_statev_offset=0,
            fp_statev_offset=80,
            allow_mock_driver=allow_mock_driver,
            name="crystal_plasticity",
            **kwargs)

    # ------------------------------------------------------------------
    # Binding construction from parsed *User Material constants
    # ------------------------------------------------------------------
    def make_binding(self, constants, name: str = "",
                     section: Optional[dict] = None) -> MaterialBinding:
        """Build a MaterialBinding from the 11 ``*User Material`` constants.

        ``constants`` is the flat list Abaqus reads after ``*User Material,
        constants=11`` (crystal type, 9 rotation entries, grain id).
        """
        c = [float(x) for x in list(constants)]
        if len(c) != self.N_PROPS:
            raise ValueError(
                "Grilli CP UMAT needs exactly %d constants "
                "(crystal type, 3x3 rotation, grain id); got %d: %r"
                % (self.N_PROPS, len(c), c))
        return MaterialBinding(material=self, constants=c,
                               n_state_vars=self.N_STATE_VARS,
                               name=name or ("grain%d" % int(round(c[self.IDX_GRAIN_ID]))),
                               section=section)

    # convenience alias
    def binding_from_user_material(self, constants, **kw) -> MaterialBinding:
        return self.make_binding(constants, **kw)

    # ------------------------------------------------------------------
    # PROPS accessors
    # ------------------------------------------------------------------
    @classmethod
    def crystal_type(cls, constants) -> int:
        return int(round(float(list(constants)[cls.IDX_CRYSTAL_TYPE])))

    @classmethod
    def crystal_type_name(cls, constants) -> str:
        return cls.CRYSTAL_TYPES.get(cls.crystal_type(constants), "unknown")

    @classmethod
    def grain_id(cls, constants):
        return float(list(constants)[cls.IDX_GRAIN_ID])

    @classmethod
    def rotation_matrix(cls, constants) -> np.ndarray:
        """3x3 crystal->sample rotation from PROPS[1:10] (row-major)."""
        c = np.asarray(list(constants), dtype=float)
        return c[cls.SLICE_ROTATION_PROPS].reshape(3, 3)

    # ------------------------------------------------------------------
    # STATEV extractors (operate on a single-IP STATEV slab of length >= 125)
    # ------------------------------------------------------------------
    @classmethod
    def extract_cauchy_stress(cls, statev) -> np.ndarray:
        """Cauchy stress (s11,s22,s33,s12,s13,s23) from SDV 48-53."""
        return np.asarray(statev, dtype=float).ravel()[cls.SLICE_STRESS].copy()

    @classmethod
    def extract_Fp(cls, statev) -> np.ndarray:
        """Plastic deformation gradient (3x3) from SDV 81-89 (row-major)."""
        return np.asarray(statev, dtype=float).ravel()[cls.SLICE_FP].reshape(3, 3)

    @classmethod
    def extract_cumulative_slip(cls, statev) -> float:
        """Cumulative plastic slip from SDV 35."""
        return float(np.asarray(statev, dtype=float).ravel()[cls.IDX_CUM_SLIP])

    @classmethod
    def extract_rotation(cls, statev) -> np.ndarray:
        """Current rotation matrix gmatinv (3x3) from SDV 1-9 (row-major)."""
        return np.asarray(statev, dtype=float).ravel()[cls.SLICE_ROTATION].reshape(3, 3)

    # ------------------------------------------------------------------
    # First-increment STATEV seeding (mirror umat.for kinc<=1 block)
    # ------------------------------------------------------------------
    def init_state(self, binding: MaterialBinding, coords=None,
                   n_ip: int = 8) -> np.ndarray:
        """Seed STATEV as the Grilli UMAT does on its first increment:
        rotation into SDV 1-9, Fp=I into SDV 81/85/89 (via the base class),
        plus SDV 54 = 0.01 (sessile SSD density).

        NOTE: this seeds only the STATEV history. The real UMAT's init block ALSO
        seeds the static /UMPS/ COMMON block (kGrainIndex, kSigma0, kFp,
        dislocation densities via kRhoTwinInit, ...); that is only reproduced by
        actually running the real UMAT at KSTEP=1, KINC=1 and cannot be set here.
        """
        state = super().init_state(binding, coords=coords, n_ip=n_ip)
        if state.shape[1] > self.IDX_SSD:
            state[:, self.IDX_SSD] = 0.01
        return state


# ---------------------------------------------------------------------------
# Offline smoke test: exercises the CP helpers, init_state seeding, and the
# python-backend evaluate path with a MOCK UMAT (NOT crystal plasticity).
#   run:  python -m residual_core.materials.crystal_plasticity_adapter
# ---------------------------------------------------------------------------
def _selftest() -> bool:
    from .umat_adapter import mock_isotropic_elastic_umat

    ok = True

    # HCP alpha-uranium example constants (Compression111.inp), crystal type 0.
    constants = [0.0, 0.89931, -0.4373, 0.0, 0.26422, 0.54336,
                 -0.79684, 0.34846, 0.71661, 0.6042, 1.0]

    cp = CrystalPlasticityAdapter(backend="python",
                                  umat_fn=mock_isotropic_elastic_umat)
    binding = cp.make_binding(constants, name="grainA")

    checks = []
    checks.append(("n_state_vars=125", binding.n_state_vars == 125))
    checks.append(("crystal type name", cp.crystal_type_name(constants) == "HCP"))
    checks.append(("grain id", cp.grain_id(constants) == 1.0))
    R = cp.rotation_matrix(constants)
    checks.append(("rotation 3x3", R.shape == (3, 3)
                   and abs(R[0, 0] - 0.89931) < 1e-12))

    # init_state seeding
    state = cp.init_state(binding, n_ip=8)
    checks.append(("init state shape", state.shape == (8, 125)))
    checks.append(("init rotation seeded",
                   np.allclose(cp.extract_rotation(state[0]), R)))
    checks.append(("init Fp = I",
                   np.allclose(cp.extract_Fp(state[0]), np.eye(3))))
    checks.append(("init SSD = 0.01", abs(state[0, cp.IDX_SSD] - 0.01) < 1e-12))

    # evaluate path with a mock (proves the CP class routes through UmatAdapter).
    # NOTE: the real CP stress needs the ifort+Abaqus UMAT; here we inject a
    # fixed-modulus elastic mock (ignoring the CP PROPS) purely to exercise the
    # STRESS/DDSDDE/STATEV plumbing offline.
    def _fixed_elastic_mock(*, dfgrd1, statev, ntens, stress, ddsdde, **kw):
        from .umat_adapter import mock_isotropic_elastic_umat
        return mock_isotropic_elastic_umat(
            dfgrd1=dfgrd1, statev=statev, ntens=ntens, stress=stress,
            ddsdde=ddsdde, props=np.array([100000.0, 0.3]))

    cp.umat_fn = _fixed_elastic_mock
    F1 = np.diag([1.0, 1.0, 1.001])
    stress, ddsdde, sv_new, diag = cp.evaluate(
        {"F0": np.eye(3), "F1": F1, "element": 1, "ip": 1},
        state[0], binding, time=(0.0, 0.0), dtime=1.0, fields=None, options=None)
    checks.append(("evaluate stress (6,)", stress.shape == (6,)))
    checks.append(("evaluate ddsdde (6,6)", ddsdde.shape == (6, 6)))
    checks.append(("evaluate statev (125,)", sv_new.shape == (125,)))
    checks.append(("Cauchy mirrored to SDV48-53",
                   np.allclose(cp.extract_cauchy_stress(sv_new), stress)))
    # uniaxial-ish stretch -> finite s33
    checks.append(("finite s33 under stretch",
                   np.isfinite(stress[2]) and abs(stress[2]) > 0.0))

    for label, passed in checks:
        ok = ok and passed
        print("  [%s] %s" % ("PASS" if passed else "FAIL", label))

    # fortran backend must refuse the mock binary (no real ifort+Abaqus build)
    cp_f = CrystalPlasticityAdapter(backend="fortran")
    refused = False
    try:
        cp_f.evaluate({"F0": np.eye(3), "F1": np.eye(3), "element": 1, "ip": 1},
                      state[0], binding, time=(0.0, 0.0), dtime=1.0,
                      fields=None, options=None)
    except RuntimeError as e:
        refused = "ifort" in str(e) and "Abaqus" in str(e)
    ok = ok and refused
    print("  [%s] fortran backend raises actionable ifort+Abaqus error"
          % ("PASS" if refused else "FAIL"))

    print("  OVERALL: %s" % ("PASS" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    import sys
    print("CrystalPlasticityAdapter offline smoke test (helpers + mock UMAT)")
    sys.exit(0 if _selftest() else 1)
