"""Run with Abaqus Python 2.7; extract final nodal fields for each step."""

import json
import sys

from odbAccess import openOdb


def vector(value):
    try:
        return list(value.dataDouble)
    except Exception:
        return list(value.data)


odb = openOdb(sys.argv[1], readOnly=True)
try:
    result = {}
    for name, step in odb.steps.items():
        frame = step.frames[-1]
        result[name] = {"time": frame.frameValue}
        for field in ("U", "RF"):
            result[name][field] = dict((str(value.nodeLabel), vector(value))
                                      for value in frame.fieldOutputs[field].values)
    with open(sys.argv[2], "w") as stream:
        json.dump(result, stream, indent=2)
finally:
    odb.close()