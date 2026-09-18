"""Direct .obj loading path (Program 2 "prescribed adapter").

Program 1 distributes ``umat_<name>_oti.obj`` (with the Fortran symbols UMAT,
UMAT_OTI_INTERNAL, UMAT_OTI_EVAL) and a completed ``umat_<name>_oti.json``
contract. ctypes cannot dlopen a relocatable ``.obj``; this module links it -- as
prescribed by the contract -- into a small loadable ``.so`` that exposes the
shared ``mat_eval_v1`` / ``mat_describe_v1`` C ABI by calling UMAT_OTI_EVAL, and
returns both the ``.so`` path and a material-package manifest so the rest of the
replay tool consumes it unchanged.

No transformation, no OTI promotion, no private source: this only *links* an
opaque object the collaborator received.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from typing import Any, Dict, Mapping, Tuple

from residual_core.runtime import (
    IncompatibleBinaryError,
    check_binary_compatibility,
    prefer_static_runtime,
    shared_library_suffix,
    static_fortran_link_flags,
)

CONTRACT_SCHEMA = "resasm_umat_oti_contract_v1"


class ObjLinkError(Exception):
    pass


def _sha(path: str) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]


# mat_eval_v1 shim that calls the object's UMAT_OTI_EVAL and returns the C-ABI
# arrays. Generated per material because the props-index -> DSIGMA_DP column map
# comes from the completed contract's parameter order.
def _shim_source(ntens: int, nprops: int, nstatev: int, nparam: int,
                 model_id: str, index_to_col) -> str:
    cases = "\n".join("         case (%d); col = %d" % (idx, col)
                      for idx, col in sorted(index_to_col.items()))
    return r"""! AUTO-GENERATED Program-2 link shim: mat_eval_v1 -> UMAT_OTI_EVAL.
module resasm_objlink_shim
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), parameter :: NT=%(NT)d, NP=%(NP)d, NS=%(NS)d, NPRM=%(NPRM)d
  character(len=*), parameter :: MODEL = "%(MODEL)s"
  type, bind(C) :: resasm_mat_desc_t
     integer(c_int) :: abi_version, ntens, nprops, nstatev, kinematics, order
  end type
contains
  function mat_describe_v1(desc_out, model_id, cap) result(rc) bind(C, name="mat_describe_v1")
    type(resasm_mat_desc_t), intent(out) :: desc_out
    character(kind=c_char), intent(out) :: model_id(*)
    integer(c_int), value :: cap
    integer(c_int) :: rc
    integer :: i, n
    desc_out%%abi_version=1; desc_out%%ntens=NT; desc_out%%nprops=NP
    desc_out%%nstatev=0; desc_out%%kinematics=0; desc_out%%order=1
    if (cap > 0) then
       n = min(len(MODEL), int(cap)-1)
       do i=1,n; model_id(i)=MODEL(i:i); end do
       model_id(n+1)=c_null_char
    end if
    rc = 0
  end function

  function mat_eval_v1(desc, props, seed, nseed, kin, dkin, st_in, dst_in, tdata, &
       stress, dstress, st_out, dst_out, ddsdde, status) result(rc) bind(C, name="mat_eval_v1")
    type(resasm_mat_desc_t), intent(in) :: desc
    real(c_double), intent(in) :: props(*)
    integer(c_int), intent(in) :: seed(*)
    integer(c_int), value :: nseed
    real(c_double), intent(in) :: kin(*)
    type(c_ptr), value :: dkin, st_in, dst_in, st_out, dst_out, status
    real(c_double), intent(in) :: tdata(*)
    real(c_double), intent(out) :: stress(*), dstress(*), ddsdde(*)
    integer(c_int) :: rc
    real*8 :: STRESS_(NT), STATEV_(max(NS,1)), DDSDDE_(NT,NT), STRAN(NT), DSTRAN(NT)
    real*8 :: TIME(2), PROPS_(NP), DSIGMA_DP(NT,NPRM), DSTATEV_DP(max(NS,1),NPRM)
    integer :: i, j, s, col, pidx
    integer(c_int), pointer :: pstat
    rc = 0
    if (desc%%abi_version/=1 .or. desc%%ntens/=NT .or. desc%%nprops/=NP) then; rc=2; call fin(); return; end if
    do s=1,nseed
      if (seed(s)<1 .or. seed(s)>NP) then; rc=4; call fin(); return; end if
    end do
    do i=1,NT; STRESS_(i)=0.d0; STRAN(i)=0.d0; DSTRAN(i)=kin(i); end do
    do i=1,max(NS,1); STATEV_(i)=0.d0; end do
    do i=1,NP; PROPS_(i)=props(i); end do
    TIME(1)=tdata(1); TIME(2)=tdata(1)
    call UMAT_OTI_EVAL(STRESS_, STATEV_, DDSDDE_, STRAN, DSTRAN, TIME, tdata(2), &
         tdata(3), tdata(4), PROPS_, NP, NT, NS, NPRM, DSIGMA_DP, DSTATEV_DP)
    do i=1,NT; stress(i)=STRESS_(i); end do
    do i=1,NT; do j=1,NT; ddsdde((i-1)*NT+j)=DDSDDE_(i,j); end do; end do
    do s=1,nseed
      pidx=seed(s); col=0
      select case (pidx)
%(CASES)s
      end select
      if (col==0) then; rc=4; call fin(); return; end if
      do i=1,NT; dstress((i-1)*nseed+s)=DSIGMA_DP(i,col); end do
    end do
    call fin()
  contains
    subroutine fin()
      if (c_associated(status)) then; call c_f_pointer(status,pstat); pstat=rc; end if
    end subroutine
  end function

  ! Regular (real) UMAT evaluation -- used ONLY by Program 2 validation to build
  ! reference primal / perturbed analyses by calling the compiled REGULAR model.
  subroutine reg_eval(props, np, dstrain, statev_in, ns, stress, ddsdde, statev_out) &
       bind(C, name="reg_eval")
    real(c_double), intent(in) :: props(*), dstrain(*), statev_in(*)
    integer(c_int), value :: np, ns
    real(c_double), intent(out) :: stress(*), ddsdde(*), statev_out(*)
    real*8 :: STRESS_(NT), STATEV_(max(ns,1)), DDSDDE_(NT,NT), STRAN(NT), DSTRAN(NT)
    real*8 :: TIME(2), PREDEF(1), DPRED(1), PROPS_(np), COORDS(3), DROT(3,3)
    real*8 :: DFGRD0(3,3), DFGRD1(3,3), DDSDDT(NT), DRPLDE(NT)
    real*8 :: SSE, SPD, SCD, RPL, DRPLDT
    character*80 CMNAME
    integer :: i, j
    SSE=0.d0; SPD=0.d0; SCD=0.d0; RPL=0.d0; DRPLDT=0.d0
    do i=1,NT; STRESS_(i)=0.d0; STRAN(i)=0.d0; DSTRAN(i)=dstrain(i)
      DDSDDT(i)=0.d0; DRPLDE(i)=0.d0; end do
    do i=1,max(ns,1); STATEV_(i)=0.d0; end do
    do i=1,ns; STATEV_(i)=statev_in(i); end do
    do i=1,np; PROPS_(i)=props(i); end do
    do i=1,3; COORDS(i)=0.d0; end do
    PREDEF(1)=0.d0; DPRED(1)=0.d0
    do i=1,3; do j=1,3; DROT(i,j)=0.d0; DFGRD0(i,j)=0.d0; DFGRD1(i,j)=0.d0; end do
      DROT(i,i)=1.d0; DFGRD0(i,i)=1.d0; DFGRD1(i,i)=1.d0; end do
    TIME(1)=0.d0; TIME(2)=0.d0
    call UMAT(STRESS_,STATEV_,DDSDDE_,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT, &
         STRAN,DSTRAN,TIME,1.d0,0.d0,0.d0,PREDEF,DPRED,CMNAME, &
         3,NT-3,NT,max(ns,1),PROPS_,np,COORDS,DROT,1.d0,1.d0,DFGRD0,DFGRD1,1,1,1,1,1,1)
    do i=1,NT; stress(i)=STRESS_(i); end do
    do i=1,NT; do j=1,NT; ddsdde((i-1)*NT+j)=DDSDDE_(i,j); end do; end do
    do i=1,ns; statev_out(i)=STATEV_(i); end do
  end subroutine
end module
""" % {"NT": ntens, "NP": nprops, "NS": nstatev, "NPRM": nparam,
       "MODEL": model_id, "CASES": cases}


def link_object(obj_path: str, contract: Mapping[str, Any], workdir: str) -> str:
    """Link the distributed .obj into a loadable .so exposing the C ABI. Returns
    the .so path. Verifies the .obj matches the completed contract's hash."""
    if contract.get("schema") != CONTRACT_SCHEMA:
        raise ObjLinkError("completed contract schema must be %r, got %r"
                           % (CONTRACT_SCHEMA, contract.get("schema")))
    declared = contract.get("object", {}).get("sha256")
    if declared and _sha(obj_path) != declared:
        raise ObjLinkError("object %r hash %r != contract object.sha256 %r"
                           % (obj_path, _sha(obj_path), declared))
    # reject a platform-incompatible object BEFORE invoking the linker, so the
    # collaborator gets a clear diagnostic instead of an opaque ld failure.
    try:
        check_binary_compatibility(obj_path, contract.get("binary"))
    except IncompatibleBinaryError as exc:
        raise ObjLinkError(str(exc)) from exc
    d = contract["dimensions"]
    ntens, nprops, nstatev, nparam = (int(d["ntens"]), int(d["nprops"]),
                                      int(d["nstatev"]), int(d["nparam"]))
    idx_to_col = {int(p["props_index"]): int(p["oti_direction"]) for p in contract["parameters"]}
    os.makedirs(workdir, exist_ok=True)
    shim = os.path.join(workdir, "objlink_shim.f90")
    open(shim, "w").write(_shim_source(ntens, nprops, nstatev, nparam,
                                       contract["model_id"], idx_to_col))
    shim_o = os.path.join(workdir, "objlink_shim.o")
    fc = os.environ.get("FC") or "gfortran"
    subprocess.check_call([fc, "-c", "-fPIC", "-ffree-form", shim, "-J" + workdir, "-o", shim_o])
    so = os.path.join(workdir, contract["model_id"] + shared_library_suffix())
    link = [fc, "-shared", shim_o, obj_path, "-o", so]
    if prefer_static_runtime():
        link = [fc, "-shared"] + static_fortran_link_flags() + [shim_o, obj_path, "-o", so]
    subprocess.check_call(link)
    return so


def package_from_contract(obj_path: str, contract: Mapping[str, Any],
                          workdir: str) -> Tuple[str, Dict[str, Any]]:
    """Link the .obj and build a resasm_material_package_v1 manifest dict that the
    existing MaterialPackage/replay path consumes. Returns (manifest_dir, manifest)."""
    so = link_object(obj_path, contract, workdir)
    d = contract["dimensions"]
    manifest = {
        "schema": "resasm_material_package_v1",
        "model_id": contract["model_id"],
        "kinematics": contract["kinematics"],
        "ntens": int(d["ntens"]), "nprops": int(d["nprops"]), "nstatev": int(d["nstatev"]),
        "parameters": [{"name": p["name"], "index": int(p["props_index"]),
                        "oti_direction": int(p["oti_direction"])}
                       for p in contract["parameters"]],
        "state_layout": [], "outputs": ["stress", "consistent_tangent", "DSIGMA_DP"],
        "sensitivity_capabilities": {"runtime_parameter_seeding": True,
                                     "maximum_directions": int(d["nparam"]), "order": 1},
        "abi": {"symbol": "mat_eval_v1", "version": 1},
        "contract_version": contract.get("contract_version"),
        "binaries": {"regular": {"hash": contract.get("regular_source_hash"), "build_id": "linked"},
                     "oti": {"path": os.path.basename(so), "hash": _sha(so), "build_id": "linked"}},
        "provenance": {"linked_from_obj": os.path.basename(obj_path)},
    }
    json.dump(manifest, open(os.path.join(workdir, "manifest.json"), "w"), indent=2)
    return workdir, manifest


def load_oti_material(obj_path: str, contract_path: str, workdir: str = None):
    """Convenience: (.obj + completed contract json) -> (MaterialPackage, base_dir)
    ready for residual_core.replay.replay_sensitivities."""
    from .package import MaterialPackage
    contract = json.load(open(contract_path))
    workdir = workdir or tempfile.mkdtemp(prefix="resasm_objlink_")
    base_dir, manifest = package_from_contract(obj_path, contract, workdir)
    return MaterialPackage(manifest, base_dir=base_dir), base_dir
