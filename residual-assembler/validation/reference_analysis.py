"""Validation-only reference analyses driven by the compiled REGULAR UMAT.

Program 2's production workflow NEVER solves a primal FE problem -- it replays a
converged analysis. This module exists ONLY for validation: it builds controlled
reference analyses (base + parameter-perturbed) by calling the *compiled regular
UMAT* (the ``reg_eval`` symbol linked from the distributed .obj) for the tangent
and stress -- never isotropic_D or any hard-coded constitutive law.

For small-strain elasticity the tangent is configuration-independent, so the
primal is a single linear solve; the point is that every stress/tangent comes
from the real regular UMAT.
"""
from __future__ import annotations

import ctypes
import os
from typing import Any, Dict, List

import numpy as np

from residual_core.core.model import Element, Model
from residual_core.core.dof_manager import DofManager
from residual_core.core.constraints import partition
from residual_core.formulations import c3d8_kernel as kern
from residual_core.replay.record import _BC
from residual_core.runtime import load_shared_library


class RegularUmat:
    """ctypes handle to the .obj's regular UMAT via the reg_eval shim symbol."""
    def __init__(self, so_path: str, ntens: int, nprops: int, nstatev: int):
        self._handle = load_shared_library(so_path)
        self.lib = self._handle.lib
        self.nt, self.np, self.ns = ntens, nprops, nstatev
        D = ctypes.c_double; IP = ctypes.POINTER(D)
        self.lib.reg_eval.argtypes = [IP, ctypes.c_int, IP, IP, ctypes.c_int, IP, IP, IP]

    def eval(self, props, dstrain, statev_in=None):
        D = ctypes.c_double
        pr = (D * self.np)(*props); ds = (D * self.nt)(*dstrain)
        si = (D * max(self.ns, 1))(*(statev_in if statev_in is not None else [0.0] * max(self.ns, 1)))
        sg = (D * self.nt)(); dd = (D * (self.nt * self.nt))(); so = (D * max(self.ns, 1))()
        self.lib.reg_eval(pr, self.np, ds, si, self.ns, sg, dd, so)
        return (np.array(sg), np.array(dd).reshape(self.nt, self.nt),
                np.array(so)[:self.ns])


# --- C3D8 benchmark meshes (B1 single element, B3 2x2x2-ish patch of 2) ------
def cube_model(drive: str, uz: float = 0.01, fz: float = 1000.0):
    nodes = {1:(0,0,0),2:(1,0,0),3:(1,1,0),4:(0,1,0),5:(0,0,1),6:(1,0,1),7:(1,1,1),8:(0,1,1),
             9:(0,0,2),10:(1,0,2),11:(1,1,2),12:(0,1,2)}
    elems = {1:[1,2,3,4,5,6,7,8], 2:[5,6,7,8,9,10,11,12]}
    m = Model(nodes={k:tuple(map(float,v)) for k,v in nodes.items()},
              elements={e:Element(e,"C3D8",list(c)) for e,c in elems.items()})
    m.element_formulation = {e:"solid_c3d8_small_strain" for e in elems}
    m.element_material = {e:"MAT" for e in elems}
    bcs = [_BC(target=n, kind="ENCASTRE") for n in (1,2,3,4)]
    if drive == "disp":
        bcs += [_BC(target=n, dof_start=3, dof_end=3, value=uz, kind="value") for n in (9,10,11,12)]
    m.boundaries = bcs
    dm = DofManager(m.nodes.keys())
    f = np.zeros(dm.ndof)
    if drive == "force":
        for n in (9,10,11,12): f[dm.node_dofs(n)[2]] += fz/4.0
    return m, dm, f, nodes, elems


def primal_solve(m, dm, reg: RegularUmat, props, fext):
    """Elastic primal solve: tangent + stress from the REGULAR UMAT (reg.eval)."""
    _, D, _ = reg.eval(props, [0.0]*reg.nt)          # elastic tangent (config-independent)
    K = np.zeros((dm.ndof, dm.ndof))
    for el in m.elements.values():
        Xe = m.coords_of(el.connectivity)
        ed = np.asarray(dm.element_dofs(el.connectivity, ("UX","UY","UZ")), int)
        K[np.ix_(ed, ed)] += kern.element_tangent(Xe=Xe, Ue=None, Dmat_ip=D, mode="small")
    fr, pi, pr = partition(m, dm); u = np.zeros(dm.ndof)
    for d, v in pr.items(): u[d] = v
    rhs = fext[fr] - (K[np.ix_(fr, pi)] @ u[pi] if pi.size else 0.0)
    u[fr] = np.linalg.solve(K[np.ix_(fr, fr)], rhs)
    return u, D


def make_record(m, dm, reg, props, fext, nodes, elems, drive, model_id, reg_hash):
    u, D = primal_solve(m, dm, reg, props, fext)
    pts = kern.ABAQUS_C3D8_GAUSS.points
    sig = {}; fint = np.zeros(dm.ndof)
    for el in m.elements.values():
        Xe = m.coords_of(el.connectivity)
        ed = np.asarray(dm.element_dofs(el.connectivity, ("UX","UY","UZ")), int)
        s = np.zeros((8, reg.nt))
        for k in range(8):
            eps = kern.b_matrix_reference(Xe, pts[k])[0] @ u[ed]
            s[k] = reg.eval(props, eps)[0]            # stress from the REGULAR UMAT
        sig[el.eid] = s
        fint[ed] += kern.element_internal_force_small_strain(Xe, s)
    react = fint - fext
    bj = [{"target":n,"kind":"ENCASTRE"} for n in (1,2,3,4)]
    cl = []
    if drive == "disp":
        bj += [{"target":n,"dof":3,"value":0.01,"kind":"value"} for n in (9,10,11,12)]
    else:
        cl = [{"node":n,"dof":3,"value":250.0} for n in (9,10,11,12)]
    return {"schema":"resasm_replay_record_v1",
            "provenance":{"model_id":model_id,"regular_hash":reg_hash,"units":"MPa,mm"},
            "kinematics":"small_strain","material":{"props":list(props)},
            "mesh":{"nodes":{str(k):list(map(float,v)) for k,v in nodes.items()},
                    "elements":{str(e):{"type":"C3D8","connectivity":list(c)} for e,c in elems.items()},
                    "node_sets":{"FIXED_FACE":[1,2,3,4],"LOADED_FACE":[9,10,11,12],
                                 "ALL_NODES":sorted(nodes.keys())},
                    "element_sets":{"DOMAIN":sorted(elems.keys())},
                    "ip_order":"abaqus_c3d8"},
            "dof":{"map":"node_major_uxuyuz"},"boundaries":bj,"loads":{"cload":cl},
            "increments":[{"time":1.0,"dtime":1.0,"u":u.tolist(),
                           "stress_ip":{str(e):s.tolist() for e,s in sig.items()},
                           "reactions":react.tolist()}]}


def response_from_record(rec, req):
    inc = rec["increments"][0]
    if req["type"] == "displacement":
        return inc["u"][3*(req["node"]-1)+req["dof"]-1]
    if req["type"] == "reaction":
        return inc["reactions"][3*(req["node"]-1)+req["dof"]-1]
    return np.array(inc["stress_ip"][str(req["element"])][req["ip"]])[req["component"]]
