"""Re-run one named recovery fixture using the producer's real verifier."""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shutil
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--umat", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    if args.work.exists():
        raise FileExistsError(f"refusing to overwrite run directory {args.work}")
    sys.path.insert(0, str(args.umat / "src"))
    sys.path.insert(0, str(args.umat / "tools"))
    from umat_oti.store import TransformStore, transform_fingerprint
    from umat_oti.abaqus.deck import generate_deck
    from transform_all import WorkItem, transform_one, put_into
    import verify_store_in_abaqus as verifier
    from export_residual_fixture import freeze

    old = json.loads(args.fixture.read_text())
    rows = [json.loads(line) for line in args.results.read_text().splitlines()
            if line.strip()]
    matches = [row for row in rows if row["source"] == old["source_id"]]
    bundled = old["source_id"] == "bundled__generic_ps/src/j2_props.f"
    if bundled:
        from residual_core.core.fixture_residual_check import BUNDLED_J2_SHA256

        source = args.umat / "UMATs/UMATs/generic_ps/j2_props.f"
        assert digest(source) == old["source_sha256"] == BUNDLED_J2_SHA256
        templates = [row for row in rows if row["source"] ==
                     "CAEAssistant-Group__UMAT-Abaqus-Isotropic-Elasticity-Isothermal-Suboutine/ISOTROPIC-ELASTICITY.for"]
        assert len(templates) == 1
        material = old["material"]
        manifest = replace(
            verifier.thaw(templates[0]["manifest"]), name="J2",
            source=Path(old["source_id"]), props=tuple(material["props"]),
            nprops=len(material["props"]), nstatv=material["nstatv"],
            material_provenance=material["provenance"],
            notes="J2 material from archived fixture; retained four-step loading; byte-identical deck required")
        prior = {"manifest": manifest.as_dict()}
    elif len(matches) != 1:
        raise ValueError("expected exactly one retained experiment for this source")
    else:
        prior = matches[0]
        source = args.cache / old["source_id"]
        assert digest(source) == old["source_sha256"] == prior["source_sha256"]
        manifest = verifier.thaw(prior["manifest"])
    assert generate_deck(manifest) == old["deck"], "loading must match the old experiment"
    assert list(manifest.props) == old["material"]["props"]
    generation = transform_fingerprint()
    recorded = json.loads((args.umat / "src/umat_oti/contract/schemas/transform_generation.json").read_text())
    assert generation == recorded["transform_fingerprint"]
    assert generation != old["transform_fingerprint"]
    args.work.mkdir(parents=True)
    if bundled:
        args.cache = args.work / "cache"
        staged_source = args.cache / old["source_id"]
        staged_source.parent.mkdir(parents=True)
        shutil.copyfile(source, staged_source)
        deck_path = args.cache / "bundled__generic_ps/decks/j2.inp"
        deck_path.parent.mkdir(parents=True)
        deck_path.write_text(old["deck"])
        source = staged_source
    item = WorkItem(old["source_id"], source, digest(source), manifest.ntens,
                    "retained author experiment and fixture tensor shape")
    transformed = transform_one(item, args.work / "transform")
    assert transformed.ok, transformed.reason
    assert transformed.metadata["compiled"], transformed.metadata
    store = TransformStore(args.work / "store")
    stored = put_into(store, item, transformed)
    jobs = []
    real_run_one = verifier.run_one

    def prefixed_run_one(**kwargs):
        logical = kwargs["job"]
        actual = "imqrf_" + logical
        kwargs["job"] = actual
        print(f"RUN {actual} in {kwargs['work_dir']}", flush=True)
        report = real_run_one(**kwargs)
        directory = kwargs["work_dir"]
        jobs.append({"job": actual, "directory": str(directory), "report": report})
        (args.work / "jobs.json").write_text(json.dumps(jobs, indent=2) + "\n")
        for path in directory.glob(actual + "*"):
            alias = directory / (logical + path.name[len(actual):])
            if not alias.exists():
                alias.symlink_to(path.name)
        return report

    verifier.run_one = prefixed_run_one
    frozen = {old["source_id"]: {
        "source_sha256": old["source_sha256"],
        "manifest": prior["manifest"], "from": str(args.results),
        "states": (prior.get("tangent") or {}).get("states", []),
    }}
    result = verifier.verify_one(
        stored, {"ntens": manifest.ntens, "form": manifest.source_form},
        {"repository": old["repository"]}, args.cache, args.work / "work",
        timeout=600, discover=False, frozen=frozen)
    result["recovery_input"] = {
        "fixture_sha256": digest(args.fixture),
        "retained_results_sha256": digest(args.results),
        "deck_byte_identical": True,
        "bundled_manifest_reconstructed": bundled,
        "numerical_evidence_reused": False,
    }
    (args.work / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result["stage"], result["reason"], flush=True)
    assert result["fingerprint"] == generation == transform_fingerprint()
    assert result["stage"] == "verified", result["reason"]
    assert all(result["evidence"].get(gate) is True for gate in (
        "abaqus_job_completed", "all_requested_outputs_present",
        "complete_history_finite", "derivatives_verified",
        "mechanically_informative", "primal_agreed")), result["evidence"]
    fixture = freeze(result, args.work / "work", len(old["original"]),
                     start=old["finite_history"]["window_starts_at_record"])
    assert fixture["deck"] == old["deck"]
    for side in ("original", "transformed"):
        before = old["finite_history"]["history_grouping"][side]
        after = fixture["finite_history"]["history_grouping"][side]
        for count in ("total_increments", "material_points_per_increment", "raw_output_records"):
            assert after[count] == before[count], (side, count, before, after)
    output = args.work / args.fixture.name
    output.write_text(json.dumps(fixture, indent=1) + "\n")
    artifacts = {str(path.relative_to(args.work)): digest(path)
                 for path in args.work.rglob("*")
                 if path.is_file() and not path.is_symlink()}
    (args.work / "hashes.json").write_text(json.dumps(artifacts, indent=2) + "\n")
    print(f"FROZEN {output} sha256={digest(output)}", flush=True)


if __name__ == "__main__":
    main()