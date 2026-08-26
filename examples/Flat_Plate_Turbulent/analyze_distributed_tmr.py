#!/usr/bin/env python3
"""Quantitative diagnostics for a distributed CBS3D SA flat-plate run.

Reads the immutable rank-local CBS partition together with one distributed VTU
snapshot and reports quantities that distinguish SA-transport errors from
post-processing errors:

* wall-resolved Cf and u_tau near the NASA station x=0.97008;
* delta_99 from the owned-node velocity profile;
* median R = nu_tilde/(kappa*u_tau*y) over 30 <= y+ <= 300;
* maximum element mu_t/mu, i.e. the value actually supplied to momentum;
* maximum projected nodal mu_t/mu written for visualisation.

No third-party Python modules are required.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

Vec = Tuple[float, float, float]
Tet = Tuple[int, int, int, int]
Face = Tuple[int, int, int]


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def dot(a: Vec, b: Vec) -> float:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def cross(a: Vec, b: Vec) -> Vec:
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a: Vec) -> float:
    return math.sqrt(dot(a, a))


def tangent(v: Vec, n: Vec) -> Vec:
    nn = norm(n)
    if not nn > 0.0:
        raise ValueError("zero normal")
    nh = (n[0]/nn, n[1]/nn, n[2]/nn)
    d = dot(v, nh)
    return (v[0]-d*nh[0], v[1]-d*nh[1], v[2]-d*nh[2])


def tet_gradients(points: Sequence[Vec], tet: Tet) -> Tuple[Vec, Vec, Vec, Vec]:
    p0, p1, p2, p3 = (points[n-1] for n in tet)
    a, b, c = sub(p1,p0), sub(p2,p0), sub(p3,p0)
    det = dot(a, cross(b,c))
    if abs(det) <= 1.0e-300:
        raise ValueError("singular TET4")
    g1 = tuple(x/det for x in cross(b,c))
    g2 = tuple(x/det for x in cross(c,a))
    g3 = tuple(x/det for x in cross(a,b))
    g0 = (-(g1[0]+g2[0]+g3[0]), -(g1[1]+g2[1]+g3[1]), -(g1[2]+g2[2]+g3[2]))
    return g0, g1, g2, g3


def velocity_gradient(points: Sequence[Vec], tet: Tet, velocity: Sequence[Vec]) -> Tuple[Vec,Vec,Vec]:
    grads = tet_gradients(points, tet)
    out: List[Vec] = []
    for comp in range(3):
        row = [0.0,0.0,0.0]
        for local,node in enumerate(tet):
            value = velocity[node-1][comp]
            for j in range(3):
                row[j] += value*grads[local][j]
        out.append((row[0],row[1],row[2]))
    return out[0],out[1],out[2]


def read_partition(rank: Path):
    bco = next(rank.glob("*.bco"))
    blines = [s.strip() for s in bco.read_text().splitlines() if s.strip() and not s.lstrip().startswith(("#","!"))]
    nflag = int(blines[0].split()[0])
    mapping = {int(p[0]):int(p[1]) for p in (line.split() for line in blines[1:1+nflag])}

    plt = next(rank.glob("*.plt"))
    with plt.open("r",encoding="utf-8") as fh:
        nelem,npoin,nboun = map(int,fh.readline().split()[:3])
        tets: List[Tet] = [(-1,-1,-1,-1)]*nelem
        for _ in range(nelem):
            p=fh.readline().split(); tets[int(p[0])-1]=tuple(map(int,p[1:5]))
        points: List[Vec] = [(0.0,0.0,0.0)]*npoin
        for _ in range(npoin):
            p=fh.readline().split(); points[int(p[0])-1]=tuple(map(float,p[1:4]))
        faces=[]
        for _ in range(nboun):
            p=fh.readline().split(); faces.append((tuple(map(int,p[:3])),int(p[3]),int(p[4])))
    return tets,points,faces,mapping


def arrays(parent: ET.Element, section: str) -> Dict[str,List[float]]:
    result: Dict[str,List[float]] = {}
    node = parent.find(section)
    if node is None:
        return result
    for a in node.findall("DataArray"):
        name=a.get("Name")
        if name and a.text:
            result[name]=[float(v) for v in a.text.split()]
    return result


def read_piece(path: Path):
    root=ET.parse(path).getroot()
    piece=root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError("VTU Piece missing: {}".format(path))
    pnode=piece.find("./Points/DataArray")
    if pnode is None or not pnode.text:
        raise ValueError("VTU points missing: {}".format(path))
    raw=[float(v) for v in pnode.text.split()]
    points=[(raw[i],raw[i+1],raw[i+2]) for i in range(0,len(raw),3)]
    pd=arrays(piece,"PointData")
    cd=arrays(piece,"CellData")
    velraw=pd.get("velocity")
    if velraw is None:
        raise ValueError("velocity missing: {}".format(path))
    velocity=[(velraw[i],velraw[i+1],velraw[i+2]) for i in range(0,len(velraw),3)]
    return points,velocity,pd,cd


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument("partition_root",type=Path)
    p.add_argument("output_dir",type=Path)
    p.add_argument("--station",type=float,default=0.97008)
    p.add_argument("--re-x1",type=float,default=5.0e6)
    p.add_argument("--u-inf",type=float,default=1.0)
    p.add_argument("--kappa",type=float,default=0.41)
    return p.parse_args()


def main():
    args=parse_args()
    ranks=sorted(p for p in args.partition_root.glob("rank_[0-9][0-9][0-9][0-9]") if p.is_dir())
    if len(ranks)!=40:
        raise SystemExit("expected 40 partition ranks, found {}".format(len(ranks)))

    rx=re.compile(r"flatplate_step_(\d{8})_rank_(\d{4})\.vtu$")
    by_step: Dict[int,Dict[int,Path]]={}
    for f in args.output_dir.glob("flatplate_step_*_rank_*.vtu"):
        m=rx.match(f.name)
        if m:
            by_step.setdefault(int(m.group(1)),{})[int(m.group(2))]=f
    good=[step for step,pieces in by_step.items() if len(pieces)==40]
    if not good:
        raise SystemExit("no complete 40-piece distributed VTU snapshot in {}".format(args.output_dir))
    step=max(good)

    nu=args.u_inf/args.re_x1
    wall=[]
    owned_nodes: Dict[int,Tuple[Vec,Vec,float,float]]={}
    max_element_ratio=0.0
    max_nodal_ratio=0.0

    for irank,rank in enumerate(ranks):
        tets,mesh_points,bfaces,mapping=read_partition(rank)
        vtu_points,velocity,pd,cd=read_piece(by_step[step][irank])
        if len(mesh_points)!=len(vtu_points):
            raise RuntimeError("mesh/VTU point-count mismatch rank {}".format(irank))

        gnode=[int(round(x)) for x in pd["global_node_id"]]
        owned=[int(round(x)) for x in pd["is_owned"]]
        nutilde=pd["nu_tilde"]
        nodalratio=pd["mu_t_over_mu"]
        for i,gid in enumerate(gnode):
            if owned[i]==1:
                owned_nodes[gid]=(vtu_points[i],velocity[i],nutilde[i],nodalratio[i])
                max_nodal_ratio=max(max_nodal_ratio,nodalratio[i])

        mu_t_e=cd["mu_t_e"]
        mu_e=cd["mu_e"]
        for ie in range(len(tets)):
            if mu_e[ie]>0.0:
                max_element_ratio=max(max_element_ratio,mu_t_e[ie]/mu_e[ie])

        for face,parent,rawbc in bfaces:
            if mapping.get(rawbc)!=530:
                continue
            pts=[mesh_points[n-1] for n in face]
            xc=sum(p[0] for p in pts)/3.0
            if xc<=0.0:
                continue
            tet=tets[parent-1]
            opp=[n for n in tet if n not in set(face)]
            if len(opp)!=1:
                raise RuntimeError("wall face/parent mismatch")
            opposite=mesh_points[opp[0]-1]
            normal=cross(sub(pts[1],pts[0]),sub(pts[2],pts[0]))
            centroid=(sum(p[0] for p in pts)/3.0,sum(p[1] for p in pts)/3.0,sum(p[2] for p in pts)/3.0)
            if dot(normal,sub(opposite,centroid))>0.0:
                normal=(-normal[0],-normal[1],-normal[2])
            area=0.5*norm(normal)
            nunit=(normal[0]/(2*area),normal[1]/(2*area),normal[2]/(2*area))
            gu=velocity_gradient(mesh_points,tet,velocity)
            traction=(nu*dot(gu[0],nunit),nu*dot(gu[1],nunit),nu*dot(gu[2],nunit))
            shear=norm(tangent(traction,nunit))
            cf=2.0*shear/(args.u_inf*args.u_inf)
            wall.append((xc,area,cf))

    if not wall or not owned_nodes:
        raise RuntimeError("incomplete wall/profile data")

    # Use a local face window determined by the nearest wall-centroid spacing.
    xs=sorted(set(round(x,12) for x,_,_ in wall))
    nearest=min(xs,key=lambda x:abs(x-args.station))
    spacings=[abs(x-nearest) for x in xs if abs(x-nearest)>1.0e-12]
    window=1.51*min(spacings) if spacings else 1.0e-8
    sample=[r for r in wall if abs(r[0]-args.station)<=window]
    if not sample:
        sample=min(wall,key=lambda r:abs(r[0]-args.station)),
    area=sum(r[1] for r in sample)
    cf=sum(r[1]*r[2] for r in sample)/area
    xcf=sum(r[1]*r[0] for r in sample)/area
    utau=args.u_inf*math.sqrt(0.5*cf)

    unique_x=sorted(set(round(v[0][0],12) for v in owned_nodes.values() if v[0][0]>=0.0))
    profile_x=min(unique_x,key=lambda x:abs(x-args.station))
    selected=[v for v in owned_nodes.values() if abs(v[0][0]-profile_x)<=2.0e-11]
    by_y: Dict[float,List[Tuple[Vec,float,float]]]={}
    for coord,vel,q,ratio in selected:
        by_y.setdefault(round(coord[1],14),[]).append((vel,q,ratio))
    profile=[]
    for y,vals in by_y.items():
        u=sum(v[0][0] for v in vals)/len(vals)
        q=sum(v[1] for v in vals)/len(vals)
        profile.append((float(y),u,q))
    profile.sort()

    delta99=float("nan")
    target=0.99*args.u_inf
    for i,(y,u,q) in enumerate(profile):
        if u>=target:
            if i==0:
                delta99=y
            else:
                y0,u0,_=profile[i-1]
                delta99=y0+(target-u0)*(y-y0)/(u-u0) if u!=u0 else y
            break

    ratios=[]
    for y,u,q in profile:
        if y<=0.0:
            continue
        yp=y*utau/nu
        if 30.0<=yp<=300.0:
            ratios.append(q/(args.kappa*utau*y))
    rmedian=float("nan")
    if ratios:
        z=sorted(ratios); m=len(z)//2
        rmedian=z[m] if len(z)%2 else 0.5*(z[m-1]+z[m])

    rex=args.re_x1*xcf
    cf_corr=0.0592*rex**(-0.2)
    result={
        "iteration":step,
        "wall_station_requested":args.station,
        "wall_station_area_weighted":xcf,
        "profile_x":profile_x,
        "cf_resolved":cf,
        "cf_power_1_5_secondary":cf_corr,
        "u_tau":utau,
        "delta99":delta99,
        "R_nu_tilde_over_kappa_utau_y_median_30_300":rmedian,
        "R_sample_count":len(ratios),
        "max_element_mu_t_over_mu_used_by_momentum":max_element_ratio,
        "max_projected_nodal_mu_t_over_mu":max_nodal_ratio,
        "owned_global_nodes":len(owned_nodes),
    }
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
