#!/usr/bin/env python3
"""G11a: independent original-U-channel first-step AFC effectiveness audit.

Compare corrected production AFC-on one-step VTU against same original G6
AFC-off one-step VTU, from identical actual step-50000 fields. This checks
the lower physical reference bound and demonstrates thermal-only changes.
One step does NOT establish convergence, mass conservation or overall accuracy.
"""

import argparse
import math
from pathlib import Path

from verify_uchannel_partition_parity_g8b import owned_fields


def insist(condition, message):
    if not condition:
        raise ValueError(message)


def max_field_delta(a, b, name):
    largest = (0.0, None)
    for gid in a:
        for x,y in zip(a[gid][name], b[gid][name]):
            d=abs(x-y)
            insist(math.isfinite(d), f"nonfinite {name} at gid {gid}")
            if d > largest[0]:
                largest=(d,gid)
    return largest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--afc-off",type=Path,required=True)
    p.add_argument("--afc-on",type=Path,required=True)
    p.add_argument("--before",type=int,default=50000)
    p.add_argument("--after",type=int,default=50001)
    p.add_argument("--reference-temperature",type=float,default=298.15)
    p.add_argument("--node",type=int,default=5707)
    p.add_argument("--bound-tol",type=float,default=1e-8)
    p.add_argument("--velocity-tol",type=float,default=1e-10)
    p.add_argument("--pressure-tol",type=float,default=1e-7)
    a=p.parse_args()
    off0=owned_fields(a.afc_off,a.before,2)
    on0=owned_fields(a.afc_on,a.before,2)
    off1=owned_fields(a.afc_off,a.after,2)
    on1=owned_fields(a.afc_on,a.after,2)
    insist(off0.keys()==on0.keys()==off1.keys()==on1.keys(),
           "global node ID sets mismatch")
    print("=== G11a ORIGINAL U-CHANNEL THERMAL AFC EFFECTIVENESS ===")
    print(f"Global owned nodes: {len(off0)}")
    for name in ("temperature","velocity","pressure","coordinates"):
        d,gid=max_field_delta(off0,on0,name)
        print(f"initial {name:12s} max difference: {d:.12e} at gid {gid}")
        insist(d<=1e-12,f"initial {name} fields differ")
    for name,tol in (("velocity",a.velocity_tol),
                     ("pressure",a.pressure_tol),
                     ("coordinates",1e-12)):
        d,gid=max_field_delta(off1,on1,name)
        print(f"final {name:12s} max difference  : {d:.12e} at gid {gid}")
        insist(d<=tol,f"unexpected flow change from thermal AFC in uncoupled case")
    toff={gid:v["temperature"][0] for gid,v in off1.items()}
    ton={gid:v["temperature"][0] for gid,v in on1.items()}
    off_min=min(toff.values())
    on_min=min(ton.values())
    off_count=sum(v < a.reference_temperature-a.bound_tol for v in toff.values())
    on_count=sum(v < a.reference_temperature-a.bound_tol for v in ton.values())
    d,gid=max_field_delta(off1,on1,"temperature")
    print(f"Tmin AFC OFF [K] : {off_min:.12f}; cold nodes={off_count}")
    print(f"Tmin AFC ON  [K] : {on_min:.12f}; cold nodes={on_count}")
    print(f"max AFC temperature effect [K]: {d:.12e} at gid {gid}")
    print(f"node {a.node}: T AFC OFF={toff[a.node]:.12f} K, "
          f"T AFC ON={ton[a.node]:.12f} K")
    print("AFC does NOT alter the pressure/velocity continuity solver.")
    insist(off_count>0,"G6 reference no longer displays original undershoot")
    insist(on_min >= a.reference_temperature-a.bound_tol,
           "AFC-on still violates reference temperature lower bound")
    print("G11A ONE-STEP THERMAL AFC REFERENCE-BOUND PASS")
    print("Caveat: global conservation and steady-state accuracy require separate tests.")


if __name__=="__main__":
    try:
        main()
    except (ValueError,OSError,IndexError) as exc:
        raise SystemExit(f"G11a AFC AUDIT ERROR: {exc}")
