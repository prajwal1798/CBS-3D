#!/usr/bin/env python3
"""CBS3D transient CHT wall-flux, HTC and Nusselt post-processor.

The audited implementation is stored as compressed payload chunks beside this
launcher to keep the GitHub write portable. Only the Python standard library is
required on Sunbird.
"""
from __future__ import print_function
import base64
import gzip
from pathlib import Path

_payload_dir = Path(__file__).resolve().parent / "_compute_cht_htc_nusselt_payload"
_parts = sorted(_payload_dir.glob("part*.txt"))
if not _parts:
    raise RuntimeError("Missing HTC/Nusselt payload files below {}".format(_payload_dir))
_payload = "".join(path.read_text(encoding="ascii").strip() for path in _parts)
_source = gzip.decompress(base64.b64decode(_payload)).decode("utf-8")
exec(compile(_source, __file__, "exec"), globals(), globals())
