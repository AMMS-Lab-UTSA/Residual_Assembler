from __future__ import print_function

import json
import sys

from odbAccess import openOdb

odb = openOdb(sys.argv[-2], readOnly=True)
instance = list(odb.rootAssembly.instances.values())[0]
step = odb.steps["LOAD"]
frames = []
for frame_index, frame in enumerate(step.frames):
    displacement = {}
    reaction = {}
    stress = {}
    for value in frame.fieldOutputs["U"].getSubset(region=instance).values:
        displacement[str(value.nodeLabel)] = [float(component) for component in value.data]
    for value in frame.fieldOutputs["RF"].getSubset(region=instance).values:
        reaction[str(value.nodeLabel)] = [float(component) for component in value.data]
    for value in frame.fieldOutputs["S"].getSubset(region=instance).values:
        stress.setdefault(str(value.elementLabel), []).append(
            [float(component) for component in value.data]
        )
    frames.append(
        {
            "frame": frame_index,
            "time": float(frame.frameValue),
            "U": displacement,
            "RF": reaction,
            "S": stress,
        }
    )
with open(sys.argv[-1], "w") as handle:
    json.dump({"frames": frames}, handle, indent=2)
odb.close()
