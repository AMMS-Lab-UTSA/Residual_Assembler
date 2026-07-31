"""Stage-2 offline replay + sensitivity (collaborator side).

This package is the collaborator's half of the two-stage, privacy-preserving
sensitivity framework:

    Stage 1 (elsewhere): run the production FE analysis with JHU's REGULAR
        material binary; save a replay record.
    Stage 2 (here):      replay that record through JHU's matched OTI-enabled
        material binary -- called only through the stable C ABI in
        ``contract/resasm_mat_abi_v1.h`` -- assemble R, K and R_,p, solve
        ``K du/dp = -R_,p``, and propagate to requested outputs dq/dp.

It never sees material source: the OTI binary is loaded through
:class:`~residual_core.replay.abi.MaterialABI` (ctypes), and the residual
math reuses the existing, tested :mod:`residual_core.core.field_sensitivity`
solver and :mod:`residual_core.formulations.c3d8_kernel`.

The concrete OTI toolchain is not required to exercise this: an isolated test
fixture (tests/fixtures/replay_elastic) can build an opaque binary implementing
the same ABI so the pipeline runs end-to-end. The material provider itself lives
in the SEPARATE Program-1 repo (UMAT_source_transformation/oti_provider); this
package never contains provider or source-transformation code.
"""

from .abi import MaterialABI, MatEvalError, ABI_VERSION, KIN_SMALL_STRAIN, KIN_FINITE_STRAIN
from .package import MaterialPackage, TwinMismatchError, PackageError
from .record import ReplayRecord, RecordError, preflight_record
from .engine import replay_sensitivities, run_request, ReplayResult
from .objlink import link_object, package_from_contract, load_oti_material, ObjLinkError
from .job import run_job, JobError

__all__ = [
    "MaterialABI", "MatEvalError", "ABI_VERSION",
    "KIN_SMALL_STRAIN", "KIN_FINITE_STRAIN",
    "MaterialPackage", "TwinMismatchError", "PackageError",
    "ReplayRecord", "RecordError", "preflight_record",
    "replay_sensitivities", "run_request", "ReplayResult",
    "link_object", "package_from_contract", "load_oti_material", "ObjLinkError",
    "run_job", "JobError",
]
