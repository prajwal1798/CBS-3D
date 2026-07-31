#!/usr/bin/env python3
"""CBS3D transient CHT wall-flux, HTC and Nusselt post-processor.

The audited numerical implementation is stored as compressed payload chunks.
Journal-only SVG restyling is implemented in a normal Python module to avoid
fragile encoded plotting payloads on older Sunbird Python installations.
"""
from __future__ import print_function

import base64
import binascii
import gzip
import sys
from pathlib import Path


def _decode_gzip_base64(text, source_path):
    """Decode whitespace-tolerant, optionally unpadded Base64."""
    compact = "".join(text.split())
    if not compact:
        raise RuntimeError("Empty compressed payload: {}".format(source_path))

    compact += "=" * ((-len(compact)) % 4)

    try:
        compressed = base64.b64decode(compact)
    except (TypeError, ValueError, binascii.Error) as exc:
        raise RuntimeError(
            "Invalid Base64 payload {}: {}".format(source_path, exc)
        )

    try:
        return gzip.decompress(compressed).decode("utf-8")
    except (IOError, OSError, EOFError, UnicodeDecodeError) as exc:
        raise RuntimeError(
            "Invalid gzip payload {}: {}".format(source_path, exc)
        )


def _run_plots_only(output_dir):
    from compute_cht_htc_nusselt_journal import regenerate_journal_plots

    regenerate_journal_plots(Path(output_dir))


if __name__ == "__main__" and len(sys.argv) == 3 \
        and sys.argv[1] == "--plots-only":
    try:
        _run_plots_only(sys.argv[2])
        sys.exit(0)
    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(1)


_tool_dir = Path(__file__).resolve().parent
_payload_dir = _tool_dir / "_compute_cht_htc_nusselt_payload"
_parts = sorted(_payload_dir.glob("part*.txt"))
if not _parts:
    raise RuntimeError(
        "Missing HTC/Nusselt payload files below {}".format(_payload_dir)
    )

_payload = "".join(
    path.read_text(encoding="ascii") for path in _parts
)
_source = _decode_gzip_base64(_payload, _payload_dir)

_original_name = globals().get("__name__", "__main__")
globals()["__name__"] = "_compute_cht_htc_nusselt_payload"
exec(compile(_source, __file__, "exec"), globals(), globals())
globals()["__name__"] = _original_name

if _original_name == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(1)
