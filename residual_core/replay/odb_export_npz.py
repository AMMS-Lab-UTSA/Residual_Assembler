# -*- coding: utf-8 -*-
"""Export every frame of a one-step ODB to the history-replay ``.npz`` layout.

Run with Abaqus Python (2.7, numpy, odbAccess):

    abaqus python odb_export_npz.py -- Analysis.odb fields.npz

Arrays (frame axis first; frame 0 is the initial state):

  time[f], increment[f]     step time and increment number of each frame
  node_labels[n], coords[n,3], elem_labels[e], conn[e,8]
  U[f,n,3], RF[f,n,3], CF[f,n,3] (CF only when it is in the ODB)
  S[f,e,8,6]                integration-point stress (11,22,33,12,13,23)
  SDV[f,e,8,nsdv]           SDV1..SDVn
  step_name, instance_name, element_type, precision

The ODB stores field output in single precision; ``precision`` records that
so the replay applies single-precision tolerances. The script refuses ODBs it
cannot export faithfully (several steps or instances, non-C3D8 elements,
missing U/RF/S at a frame) instead of writing a partial file.
"""
import sys

import numpy as np
from odbAccess import openOdb


def fail(message):
    sys.stderr.write("odb_export_npz: %s\n" % message)
    sys.exit(2)


def nodal(field, index, count, frame_number, name):
    out = np.zeros((count, 3))
    seen = np.zeros(count, dtype=bool)
    for block in field.bulkDataBlocks:
        rows = np.array([index[label] for label in block.nodeLabels])
        data = np.asarray(block.data, dtype=float)
        out[rows, :data.shape[1]] = data
        seen[rows] = True
    if not seen.all():
        fail("frame %d: %s is missing at %d nodes" % (frame_number, name, int((~seen).sum())))
    return out


def points(field, index, count, width, frame_number, name):
    out = np.zeros((count, 8, width))
    seen = np.zeros((count, 8), dtype=bool)
    for block in field.bulkDataBlocks:
        rows = np.array([index[label] for label in block.elementLabels])
        ips = np.asarray(block.integrationPoints, dtype=int) - 1
        data = np.asarray(block.data, dtype=float).reshape(len(rows), -1)
        out[rows, ips, :] = data[:, :width]
        seen[rows, ips] = True
    if not seen.all():
        fail("frame %d: %s is missing at %d integration points" % (frame_number, name, int((~seen).sum())))
    return out


def main(argv):
    args = argv[argv.index('--') + 1:] if '--' in argv else argv[1:]
    if len(args) != 2:
        fail("usage: abaqus python odb_export_npz.py -- Analysis.odb fields.npz")
    odb = openOdb(args[0], readOnly=True)
    if len(odb.steps.keys()) != 1:
        fail("exactly one step is supported (found %d)" % len(odb.steps.keys()))
    instances = [name for name in odb.rootAssembly.instances.keys() if name != 'ASSEMBLY']
    if len(instances) != 1:
        fail("exactly one part instance is supported (found %s)" % instances)
    instance = odb.rootAssembly.instances[instances[0]]
    nodes = sorted(instance.nodes, key=lambda n: n.label)
    node_labels = np.array([n.label for n in nodes])
    coords = np.array([n.coordinates for n in nodes], dtype=float)
    nidx = dict((label, i) for i, label in enumerate(node_labels))
    elements = sorted(instance.elements, key=lambda e: e.label)
    types = sorted(set(e.type for e in elements))
    if types != ['C3D8']:
        fail("only C3D8 elements are supported (found %s)" % types)
    elem_labels = np.array([e.label for e in elements])
    conn = np.array([e.connectivity for e in elements])
    eidx = dict((label, i) for i, label in enumerate(elem_labels))
    step_name = odb.steps.keys()[0]
    frames = odb.steps[step_name].frames
    last = frames[-1].fieldOutputs
    nsdv = len([key for key in last.keys() if key.startswith('SDV') and key[3:].isdigit()])
    nf, nn, ne = len(frames), len(nodes), len(elements)
    arrays = dict(U=np.zeros((nf, nn, 3)), RF=np.zeros((nf, nn, 3)), S=np.zeros((nf, ne, 8, 6)),
                  SDV=np.zeros((nf, ne, 8, nsdv)))
    has_cf = 'CF' in last.keys()
    if has_cf:
        arrays['CF'] = np.zeros((nf, nn, 3))
    time = np.zeros(nf)
    increment = np.zeros(nf, dtype=int)
    for f, frame in enumerate(frames):
        time[f] = frame.frameValue
        increment[f] = frame.incrementNumber
        fo = frame.fieldOutputs
        for key in ('U', 'RF', 'S'):
            if key not in fo.keys():
                fail("frame %d has no %s field output" % (f, key))
        arrays['U'][f] = nodal(fo['U'], nidx, nn, f, 'U')
        arrays['RF'][f] = nodal(fo['RF'], nidx, nn, f, 'RF')
        if has_cf and 'CF' in fo.keys():
            for block in fo['CF'].bulkDataBlocks:
                rows = np.array([nidx[label] for label in block.nodeLabels])
                arrays['CF'][f][rows] = np.asarray(block.data, dtype=float)
        arrays['S'][f] = points(fo['S'], eidx, ne, 6, f, 'S')
        for s in range(nsdv):
            key = 'SDV%d' % (s + 1)
            if key not in fo.keys():
                if f == 0:
                    continue          # Abaqus writes no SDV at the virgin frame: zero initial state
                fail("frame %d has no %s" % (f, key))
            arrays['SDV'][f][:, :, s] = points(fo[key], eidx, ne, 1, f, key)[:, :, 0]
    precision = 'float32'
    try:
        if str(frames[-1].fieldOutputs['U'].values[0].precision) == 'DOUBLE_PRECISION':
            precision = 'float64'
    except (AttributeError, IndexError):
        pass
    np.savez_compressed(args[1], time=time, increment=increment, node_labels=node_labels,
                        coords=coords, elem_labels=elem_labels, conn=conn, step_name=step_name,
                        instance_name=instances[0], element_type='C3D8', precision=precision, **arrays)
    odb.close()
    print("exported %d frames, %d nodes, %d elements, %d SDV, CF=%s, precision=%s -> %s"
          % (nf, nn, ne, nsdv, has_cf, precision, args[1]))


if __name__ == '__main__':
    main(sys.argv)
