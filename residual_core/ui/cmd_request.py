"""Presentation-matching collaborator CLI."""

import sys

from ..replay.presentation import RequestFailure, run_request


def register(subparsers):
    parser = subparsers.add_parser("request", help="sensitivities from Analysis.inp, Analysis.odb, OTI_UMAT.obj and sensitivity_request.json")
    for name in ("model", "odb", "material", "request", "out"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--mapping", help="advanced: generated Mapping.json; default discover beside object")
    parser.add_argument("--abaqus", default="abaqus", help="advanced: licensed Abaqus executable for ODB extraction")
    parser.add_argument("--validate", action="store_true", help="explicit independent ORIGINAL finite-difference validation")
    # Models inside the bounded presentation scope keep this module's engine;
    # anything else (prescribed displacements, many increments, any provider)
    # goes to the history replay engine. route_request decides and says which.
    from .cmd_history import route_request
    parser.set_defaults(func=route_request)


def run(args):
    try:
        result = run_request(**{key: getattr(args, key) for key in (
            "model", "odb", "material", "request", "out", "mapping", "abaqus", "validate")})
        print("request executed: %d scalar results; verified=%s" % (len(result["results"]), result["metadata"]["verified"]))
        return 0
    except RequestFailure as error:
        print("request failed: %s" % error, file=sys.stderr)
        return 2
    except (ValueError, OSError):
        print("request failed: Category: output_directory\n"
              "Action: use a new or empty writable output directory; refusing stale public results.\n"
              "Private diagnostics: unavailable (output directory could not be prepared or written).", file=sys.stderr)
        return 2