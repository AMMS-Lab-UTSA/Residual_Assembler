# -*- coding: utf-8 -*-
"""Export every frame of a cantilever ODB to a compact .npz (run with `abaqus python`).

    abaqus python export_odb.py -- claude_j2_nominal.odb claude_j2_nominal_fields.npz

Arrays (frame axis first; frame 0 is the initial state):
  time[f]            step time of the frame
  node_labels[n], coords[n,3]
  elem_labels[e], conn[e,8]
  U[f,n,3], RF[f,n,3]
  S[f,e,8,6]         integration-point stress, Abaqus order 11,22,33,12,13,23
  E[f,e,8,6]         integration-point total strain (engineering shear)
  SDV[f,e,8,nsdv]
  tipmid_U2_history[k] (history output, as stored)
Field output in an ODB is single precision; that limits any finite difference
taken from these files and is recorded in the manifest written alongside.
"""
import sys
import numpy as np
from odbAccess import openOdb

args = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else sys.argv[1:]
odb_path, out_path = args[0], args[1]
odb = openOdb(odb_path, readOnly=True)
inst = odb.rootAssembly.instances[odb.rootAssembly.instances.keys()[0]]
nodes = sorted(inst.nodes, key=lambda n: n.label)
node_labels = np.array([n.label for n in nodes])
coords = np.array([n.coordinates for n in nodes], dtype=float)
nidx = dict((lab, i) for i, lab in enumerate(node_labels))
elems = sorted(inst.elements, key=lambda e: e.label)
elem_labels = np.array([e.label for e in elems])
conn = np.array([e.connectivity for e in elems])
eidx = dict((lab, i) for i, lab in enumerate(elem_labels))
step = odb.steps[odb.steps.keys()[-1]]
frames = step.frames
nf, nn, ne = len(frames), len(nodes), len(elems)
nsdv = len([k for k in frames[-1].fieldOutputs.keys() if k.startswith('SDV')])
U = np.zeros((nf, nn, 3)); RF = np.zeros((nf, nn, 3))
S = np.zeros((nf, ne, 8, 6)); E = np.zeros((nf, ne, 8, 6)); SDV = np.zeros((nf, ne, 8, nsdv))
time = np.zeros(nf)
for f, fr in enumerate(frames):
    time[f] = fr.frameValue
    fo = fr.fieldOutputs
    for key, arr in (('U', U), ('RF', RF)):
        if key in fo.keys():
            for v in fo[key].values:
                arr[f, nidx[v.nodeLabel], :] = v.data
    for key, arr in (('S', S), ('E', E)):
        if key in fo.keys():
            for v in fo[key].values:
                arr[f, eidx[v.elementLabel], v.integrationPoint - 1, :] = v.data
    for s in range(nsdv):
        key = 'SDV%d' % (s + 1)
        if key in fo.keys():
            for v in fo[key].values:
                SDV[f, eidx[v.elementLabel], v.integrationPoint - 1, s] = v.data
hist = []
for rname, reg in step.historyRegions.items():
    if 'U2' in reg.historyOutputs.keys():
        hist = [val for _, val in reg.historyOutputs['U2'].data]
np.savez_compressed(out_path, time=time, node_labels=node_labels, coords=coords,
                    elem_labels=elem_labels, conn=conn, U=U, RF=RF, S=S, E=E, SDV=SDV,
                    tipmid_U2_history=np.array(hist))
odb.close()
print('exported %d frames, %d nodes, %d elements, %d SDV -> %s' % (nf, nn, ne, nsdv, out_path))
