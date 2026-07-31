#!/usr/bin/env python3
"""CBS3D transient CHT wall-flux, HTC and Nusselt post-processor.

The audited numerical implementation is stored as compressed payload chunks.
A separate compressed plotting layer supplies publication-quality SVG figures
without altering any numerical result. Only the Python standard library is
required on Sunbird.
"""
from __future__ import print_function

import base64
import gzip
import sys
from pathlib import Path


_tool_dir = Path(__file__).resolve().parent
_payload_dir = _tool_dir / "_compute_cht_htc_nusselt_payload"
_parts = sorted(_payload_dir.glob("part*.txt"))
if not _parts:
    raise RuntimeError(
        "Missing HTC/Nusselt payload files below {}".format(_payload_dir)
    )

_payload = "".join(
    path.read_text(encoding="ascii").strip() for path in _parts
)
_source = gzip.decompress(base64.b64decode(_payload)).decode("utf-8")

_plot_payload_path = (
    _tool_dir / "_compute_cht_htc_nusselt_journal_plot.txt"
)
if not _plot_payload_path.is_file():
    raise RuntimeError(
        "Missing publication plotting payload: {}".format(
            _plot_payload_path
        )
    )
_plot_payload = _plot_payload_path.read_text(
    encoding="ascii"
).strip()
_plot_source = gzip.decompress(
    base64.b64decode(_plot_payload)
).decode("utf-8")

_original_name = globals().get("__name__", "__main__")
globals()["__name__"] = "_compute_cht_htc_nusselt_payload"
exec(compile(_source, __file__, "exec"), globals(), globals())
exec(
    compile(
        _plot_source,
        str(_plot_payload_path),
        "exec",
    ),
    globals(),
    globals(),
)
globals()["__name__"] = _original_name

if _original_name == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(1)
