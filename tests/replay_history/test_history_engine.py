"""The history engine against independent references.

References: exact elastic scaling laws; whole-model central finite
differences of the ORIGINAL UMAT re-equilibrated in Python (a step ladder
with a plateau); the residual's own finite-difference Jacobian; a replay of a
Python-equilibrated history (a recording without single-precision loss).
"""
import numpy as np
import pytest

from residual_core.replay.history import HistoryEngine, ReplayMismatch, run_history, von_mises
from residual_core.replay.history_inputs import RecordedFields, UnsupportedFeature, read_history_model
from residual_core.replay.history_verify import summarize_fd, tangent_check, whole_model_fd

from conftest import beam_model

ELASTIC = (200000.0, 0.3, 1e9, 2000.0)          # yield stress out of reach: linear elastic J2


def test_elastic_displacement_control_scaling(provider_factory, tmp_path, times):
    """Prescribed displacements only: u is E-independent, S and RF scale with E."""
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, props=ELASTIC, steps=3)
    engine = HistoryEngine(model, material)
    result = run_history(engine, times=times(3), rtol=1e-13)
    E = ELASTIC[0]
    for record in result.increments:
        assert np.abs(record.u).max() > 0
        assert np.abs(record.du[:, 0]).max() * E <= 1e-12 * np.abs(record.u).max()
        np.testing.assert_allclose(record.dstress[..., 0] * E, record.stress,
                                   atol=1e-12 * np.abs(record.stress).max())
        reaction = record.reaction[engine.constrained]
        np.testing.assert_allclose(record.dreaction[engine.constrained, 0] * E, reaction,
                                   atol=1e-12 * np.abs(reaction).max())
        assert np.abs(record.dstress[..., 2:]).max() == 0.0      # SIGY0, H: no plastic branch


def test_elastic_load_control_scaling(provider_factory, tmp_path, times):
    """Concentrated load only: du/dE = -u/E and the stresses do not depend on E."""
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, props=ELASTIC, steps=3, push=0.0, cload=-400.0)
    engine = HistoryEngine(model, material)
    result = run_history(engine, times=times(3), rtol=1e-13)
    E = ELASTIC[0]
    for record in result.increments:
        np.testing.assert_allclose(record.du[:, 0] * E, -record.u, atol=1e-12 * np.abs(record.u).max())
        assert np.abs(record.dstress[..., 0]).max() * E <= 1e-11 * np.abs(record.stress).max()
    assert np.abs(model.load_end).sum() == 400.0


def _fd_summary(engine, steps, ladder, rtol=1e-13):
    result = whole_model_fd(engine, np.linspace(0, 1, steps + 1), steps=ladder, rtol=rtol)
    return result, summarize_fd(result)


def _assert_fd(summary, tolerance, expected_nonzero=None):
    for name, rows in summary.items():
        for key, row in rows.items():
            assert row["zero_reference_max_weighted_error"] <= 1e-7, (name, key, row)
            if row["nonzero_increments"]:
                assert row["max_error"] <= tolerance, (name, key, row)
                assert row["max_error"] <= 3 * row["max_spread"] + 1e-9, (name, key, row)
    if expected_nonzero:
        for name, count in expected_nonzero.items():
            assert summary[name]["U"]["nonzero_increments"] == count, name


def test_j2_whole_model_fd_of_original_umat(provider_factory, tmp_path):
    """12x4x2 J2 beam into plasticity: every derivative within 1e-6 of FD (plateau)."""
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, n=(12, 4, 2), push=0.08, steps=10)
    engine = HistoryEngine(model, material)
    result, summary = _fd_summary(engine, 10, (1e-3, 3e-4, 1e-4, 3e-5, 1e-5))
    reference = result["reference_result"]
    plastic = [int((r.state > 0).sum()) for r in reference.increments]
    assert plastic[3] == 0 and plastic[4] > 0 and plastic[-1] > 100    # elastic, then yielding
    _assert_fd(summary, 1e-6, expected_nonzero={"SIGY0": 6, "H": 6, "nu": 10})


def test_fcc_whole_model_fd_of_original_umat(provider_factory, tmp_path):
    """Crystal plasticity, 10 parameters in provider order, rate-dependent (dt matters)."""
    _, _, material = provider_factory("m6_fcc")
    props = (168000.0, 121000.0, 75000.0, 13.0, 55.0, 800.0, 2.0, 1.4, 0.001, 0.05)
    _, model = beam_model(tmp_path, n=(6, 2, 1), push=0.012, steps=5, props=props, depvar=12)
    engine = HistoryEngine(model, material)
    result, summary = _fd_summary(engine, 5, (1e-3, 1e-4, 1e-5, 1e-6))
    hardened = result["reference_result"].increments[-1].state.max()
    assert hardened > 1e-2                                     # slip hardening is active
    _assert_fd(summary, 1e-6)
    assert all(summary[name]["MISES"]["nonzero_increments"] == 5 for name in ("C11", "C12", "C44"))
    assert all(summary[name]["MISES"]["nonzero_increments"] >= 3 for name in material.params)


def test_state_that_reads_total_strain_whole_model_fd(provider_factory, tmp_path):
    """The damage fixture uses STRAN and a history maximum: the STRAN seed is exercised.

    With A < 0 its law sigma = E K0 exp(-A (eps - K0)) hardens (A > 0 softens,
    which would make the displacement-driven beam problem ill-posed).
    """
    _, _, material = provider_factory("damage")
    props = (70000.0, 0.25, 2e-4, -400.0, 0.05)
    _, model = beam_model(tmp_path, n=(6, 2, 1), push=0.006, steps=5, props=props)
    engine = HistoryEngine(model, material)
    result, summary = _fd_summary(engine, 5, (1e-3, 1e-4, 1e-5))
    assert result["reference_result"].increments[-1].state.max() > 2e-4        # damaged
    _assert_fd(summary, 1e-6)


def test_stiffness_is_the_jacobian_of_the_original_residual(provider_factory, tmp_path, times):
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, n=(4, 2, 1), push=0.04, steps=2)
    engine = HistoryEngine(model, material)
    result = run_history(engine, times=times(2), rtol=1e-13)
    last = result.increments[-1]
    assert (last.state > 0).any()
    previous = result.increments[0]
    npts = engine.ne * 8
    strain_prev = engine.flat(engine.strain(previous.u))
    common = dict(time=np.array([0.5, 0.5]), dtime=0.5, coords=engine.ip_coords.reshape(-1, 3),
                  celent=engine.celent, noel=engine.noel, npt=engine.npt)
    zero = np.zeros((npts, 6, 4))

    def residual(u):
        dstran = engine.flat(engine.strain(u)) - strain_prev
        out = material.regular(model.props, previous.stress.reshape(npts, 6), previous.state.reshape(npts, 1),
                               strain_prev, dstran, **common)
        return engine.force(engine.unflat(out["stress"]))

    dstran = engine.flat(engine.strain(last.u)) - strain_prev
    tangent = material.total(model.props, previous.stress.reshape(npts, 6), previous.state.reshape(npts, 1),
                             strain_prev, dstran, dstress_in=zero, dstate_in=np.zeros((npts, 1, 4)),
                             stran_dp=zero, dstran_dp=zero, **common)["ddsdde"]
    stiffness = engine.stiffness(engine.unflat(tangent)).toarray()
    for dof in engine.free[:: max(1, len(engine.free) // 12)]:
        columns = []
        for h in (2e-9, 1e-9):
            shift = np.zeros(engine.ndof)
            shift[dof] = h
            columns.append((residual(last.u + shift) - residual(last.u - shift)) / (2 * h))
        scale = np.abs(columns[1]).max()
        assert np.abs(columns[0] - columns[1]).max() / scale < 1e-5
        assert np.abs(stiffness[:, dof] - columns[1]).max() / scale < 1e-5


def _record_from(result, model, engine):
    count = len(result.increments) + 1
    U = np.zeros((count, engine.ndof))
    RF = np.zeros((count, engine.ndof))
    S = np.zeros((count, engine.ne, 8, 6))
    SDV = np.zeros((count, engine.ne, 8, engine.material.nstatev))
    for n, record in enumerate(result.increments, 1):
        U[n], S[n], SDV[n] = record.u, record.stress, record.state
        RF[n, engine.constrained] = record.reaction[engine.constrained]
    nn = len(model.node_ids)
    return RecordedFields(time=result.times, U=U.reshape(count, nn, 3), RF=RF.reshape(count, nn, 3),
                          S=S, SDV=SDV, source="python", precision="float64")


def test_replay_of_an_equilibrated_record_reproduces_the_solve(provider_factory, tmp_path, times):
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, n=(8, 2, 2), push=0.06, steps=6)
    engine = HistoryEngine(model, material)
    solved = run_history(engine, times=times(6), rtol=1e-13)
    replayed = run_history(engine, fields=_record_from(solved, model, engine))
    for a, b in zip(solved.increments, replayed.increments):
        np.testing.assert_allclose(b.du, a.du, rtol=0, atol=1e-10 * np.abs(a.du).max())
        np.testing.assert_allclose(b.dstress, a.dstress, rtol=0, atol=1e-10 * np.abs(a.dstress).max())
        assert b.parity["stress_max_abs"] <= 1e-9 * np.abs(a.stress).max()
    assert replayed.mode == "replay"
    check = tangent_check(engine, replayed, [1, 6], sample=8)
    assert check["max_relative_error"] < 1e-5 and check["fd_plateau_spread"] < 1e-5


def test_replay_mismatches_are_failures(provider_factory, tmp_path, times):
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, n=(6, 2, 2), push=0.05, steps=4)
    engine = HistoryEngine(model, material)
    solved = run_history(engine, times=times(4), rtol=1e-13)
    stress_off = _record_from(solved, model, engine)
    stress_off.S[3, 2, 5, 0] += 1e-3 * np.abs(stress_off.S[3]).max()
    with pytest.raises(ReplayMismatch, match="increment 3 .*stress"):
        run_history(engine, fields=stress_off)
    moved = _record_from(solved, model, engine)
    moved.U[2, 10, 1] += 1e-4 * np.abs(moved.U[2]).max()          # a free node off equilibrium
    with pytest.raises(ReplayMismatch, match="increment 2"):
        run_history(engine, fields=moved)
    other = _record_from(solved, model, engine)
    stiffer = model.props.copy()
    stiffer[0] *= 1.01
    with pytest.raises(ReplayMismatch):
        run_history(engine, fields=other, props=stiffer)


def test_reequilibration_polishes_to_double_precision(provider_factory, tmp_path, times):
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, n=(8, 2, 2), push=0.06, steps=5)
    engine = HistoryEngine(model, material)
    solved = run_history(engine, times=times(5), rtol=1e-13)
    record = _record_from(solved, model, engine)
    record.U = record.U.astype(np.float32).astype(np.float64)     # what an ODB stores
    record.U.reshape(len(times(5)), -1)[:, model.prescribed_dofs] = np.array(
        [model.prescribed_at(t) for t in times(5)])
    rough = run_history(engine, fields=record)
    polished = run_history(engine, fields=record, reequilibrate=True)
    for a, b, c in zip(solved.increments, rough.increments, polished.increments):
        assert c.residual_free_max <= 1e-11 * c.residual_scale
        assert b.residual_free_max > c.residual_free_max
        np.testing.assert_allclose(c.du, a.du, rtol=0, atol=1e-9 * np.abs(a.du).max())


@pytest.mark.parametrize("change, message", [
    (lambda t: t.replace("nlgeom=NO", "nlgeom=YES"), "NLGEOM"),
    (lambda t: t.replace("type=C3D8,", "type=C3D8R,"), "C3D8R"),
    (lambda t: t.replace("*Boundary\nTIP", "*Amplitude, name=A1\n0., 0., 1., 1.\n*Boundary\nTIP"), "Amplitude"),
    (lambda t: t.replace("*Boundary\nTIP", "*Boundary, op=NEW\nTIP"), "OP=NEW"),
    (lambda t: t.replace("*Boundary\nROOT, ENCASTRE", "*Boundary\nROOT, 1, 1, 0.01"), "model-data"),
    (lambda t: t.replace("*End Step", "*Dsload\nALL, P, 1.0\n*End Step"), "Dsload"),
])
def test_unsupported_deck_features_are_named(tmp_path, change, message):
    with pytest.raises(UnsupportedFeature, match=message):
        beam_model(tmp_path, extra=change)


def test_depvar_and_nprops_must_match_the_provider(provider_factory, tmp_path):
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, depvar=2)
    with pytest.raises(UnsupportedFeature, match="Depvar"):
        HistoryEngine(model, material)
    _, model = beam_model(tmp_path, props=(1.0, 0.3, 250.0, 2000.0, 7.0), name="five.inp")
    with pytest.raises(UnsupportedFeature, match="NPROPS"):
        HistoryEngine(model, material)
