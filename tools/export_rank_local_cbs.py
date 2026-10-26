#!/usr/bin/env python3
# ================================
# FILE      : export_rank_local_cbs.py
# PURPOSE   : SPLIT ONE VALIDATED SERIAL CBS3D CASE INTO RANK-LOCAL CASES
#             FOR cbs3d_parallel (MPI DOMAIN DECOMPOSITION).
#
# INPUT     : <case>.plt  <case>.bco  <case>.material   (SERIAL CASE)
#             <case>.par  <case>.matprop                (COPIED TO EVERY RANK)
#             ONE PARTITION ID PER TETRAHEDRON, FROM ONE OF:
#               (A) GMSH METIS, RUN BY THIS SCRIPT        (DEFAULT)
#               (B) --partition-file  ONE INTEGER PER LINE, 0-BASED RANK
#
# OUTPUT    : <out>/rank_XXXX/<name>_rank_XXXX.{plt,bco,material,par,matprop,mpi}
#             <out>/partition_manifest.json
#             <out>/partition_counts.csv
#
# RULES     : (THE SAME AS THE ORIGINAL TOWER EXPORTER, CHECKED BYTE FOR BYTE
#              AGAINST cases/flat_plate_p8)
#             1. A NODE IS OWNED BY THE LOWEST RANK THAT HAS AN ELEMENT ON IT.
#             2. LOCAL NODES  = OWNED NODES (ASCENDING GLOBAL ID), THEN GHOST
#                               NODES (ASCENDING GLOBAL ID).
#             3. LOCAL ELEMS  = OWNED TETRAHEDRA, ASCENDING GLOBAL ID.
#                               NO GHOST ELEMENTS.
#             4. BOUNDARY     = SERIAL BOUNDARY FACES WHOSE PARENT ELEMENT IS
#                               ON THIS RANK, IN SERIAL ORDER.
#             5. SEND TO q    = MY OWNED NODES THAT ARE GHOSTS ON q.
#                RECV FROM q  = MY GHOST NODES OWNED BY q.
#                BOTH LISTS IN ASCENDING LOCAL ID. A NEIGHBOUR IS LISTED ONLY
#                WHEN ITS SEND OR RECV LIST IS NON-EMPTY.
#             6. GLOBAL IDS IN THE .mpi FILE ARE THE SERIAL CASE IDS (1-BASED).
#
# USAGE     : python3 tools/export_rank_local_cbs.py SERIAL_BASE --parts N \
#                 --out-dir DIR [--case-name NAME] [--partition-file F]
#                 [--gmsh gmsh] [--overwrite]
#
# STATUS    : PYTHON 3.6 COMPATIBLE (SUNBIRD python3). NO NUMPY. NO GMSH
#             PYTHON MODULE. THE GMSH EXECUTABLE IS ONLY NEEDED FOR (A).
# ================================

from __future__ import print_function

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile


# ONE REAL NUMBER. THE .plt WRITER LEAVES NO SPACE BEFORE A NEGATIVE NUMBER,
# SO COORDINATES ARE FOUND BY PATTERN, NOT BY SPLITTING ON SPACES.
REAL = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


# ================================
# READ THE SERIAL CASE
# ================================
def read_serial_case(base):
    # .plt : HEADER, ELEMENTS, NODES, BOUNDARY FACES
    with open(base + ".plt", "r") as f:
        header = f.readline().split()
        nelem = int(header[0])
        npoin = int(header[1])
        nboun = int(header[2])

        elements = [None] * (nelem + 1)
        for ie in range(1, nelem + 1):
            fields = f.readline().split()
            if int(fields[0]) != ie:
                raise SystemExit("ERROR: .plt ELEMENT IDS ARE NOT 1..nelem AT %d" % ie)
            elements[ie] = (int(fields[1]), int(fields[2]),
                            int(fields[3]), int(fields[4]))

        # NODE ID IS THE FIRST 10 CHARACTERS. A NEGATIVE x CAN TOUCH THE ID,
        # SO THE COORDINATE TEXT IS KEPT EXACTLY AS WRITTEN.
        coord_text = [None] * (npoin + 1)
        for ip in range(1, npoin + 1):
            line = f.readline().rstrip("\n")
            if int(line[:10]) != ip:
                raise SystemExit("ERROR: .plt NODE IDS ARE NOT 1..npoin AT %d" % ip)
            coord_text[ip] = line[10:]

        faces = []
        for ib in range(nboun):
            fields = f.readline().split()
            faces.append((int(fields[0]), int(fields[1]), int(fields[2]),
                          int(fields[3]), int(fields[4])))

    # .bco : FIRST LINE nflag mass_flow_flag
    with open(base + ".bco", "r") as f:
        mass_flow_flag = int(f.readline().split()[1])

    # .material : id n1 n2 n3 n4 material_id
    material = [0] * (nelem + 1)
    with open(base + ".material", "r") as f:
        for ie in range(1, nelem + 1):
            fields = f.readline().split()
            if int(fields[0]) != ie:
                raise SystemExit("ERROR: .material IDS ARE NOT 1..nelem AT %d" % ie)
            material[ie] = int(fields[5])

    return {
        "nelem": nelem, "npoin": npoin, "nboun": nboun,
        "elements": elements, "coord_text": coord_text, "faces": faces,
        "mass_flow_flag": mass_flow_flag, "material": material,
    }


# ================================
# PARTITION WITH GMSH (METIS)
# ================================
# WRITES A TET-ONLY MSH2 FILE WHOSE NODE AND ELEMENT TAGS ARE THE CBS IDS,
# RUNS gmsh -part, AND READS THE PARTITION TAG BACK BY ELEMENT TAG.
def partition_with_gmsh(case, parts, gmsh_exe, work_dir):
    source = os.path.join(work_dir, "cbs_tets.msh")
    result = os.path.join(work_dir, "cbs_tets_partitioned.msh")

    with open(source, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n%d\n" % case["npoin"])
        for ip in range(1, case["npoin"] + 1):
            xyz = REAL.findall(case["coord_text"][ip])
            if len(xyz) != 3:
                raise SystemExit("ERROR: CANNOT READ COORDINATES OF NODE %d" % ip)
            f.write("%d %s %s %s\n" % (ip, xyz[0], xyz[1], xyz[2]))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % case["nelem"])
        for ie in range(1, case["nelem"] + 1):
            n = case["elements"][ie]
            f.write("%d 4 2 1 1 %d %d %d %d\n" % (ie, n[0], n[1], n[2], n[3]))
        f.write("$EndElements\n")

    options = ["-part", str(parts),
               "-format", "msh2",
               "-setnumber", "Mesh.PartitionOldStyleMsh2", "1",
               "-setnumber", "Mesh.PartitionCreateGhostCells", "0",
               "-setnumber", "Mesh.PartitionSplitMeshFiles", "0",
               "-preserve_numbering_msh2", "-save_all",
               "-o", result]

    # GMSH 4.12+ NEEDS -save TO WRITE A RE-LOADED MESH; THE OLD SUNBIRD
    # MODULE (4.4) WROTE IT WITHOUT. TRY BOTH, KEEP THE FIRST THAT WORKS.
    log_path = os.path.join(work_dir, "gmsh_partition.log")
    for extra in (["-save"], []):
        if os.path.isfile(result):
            os.remove(result)
        command = [gmsh_exe, source] + extra + options
        print("RUNNING: " + " ".join(command))
        with open(log_path, "a") as log:
            subprocess.check_call(command, stdout=log, stderr=subprocess.STDOUT)
        if os.path.isfile(result) and os.path.getsize(result) > 0:
            return read_msh2_partition(result, case["nelem"], parts)

    raise SystemExit("ERROR: GMSH WROTE NO PARTITIONED MESH; SEE " + log_path)


# ================================
# READ PARTITION IDS FROM AN OLD-STYLE PARTITIONED MSH2
# ================================
# ELEMENT LINE: id type ntags phys geom npart part1 [ghost ...] nodes...
# GMSH PARTITION IDS ARE 1..N; MPI RANK = PARTITION ID - 1.
def read_msh2_partition(path, nelem, parts):
    rank_of = [-1] * (nelem + 1)
    with open(path, "r") as f:
        for line in f:
            if line.strip() == "$Elements":
                break
        count = int(f.readline())
        for i in range(count):
            fields = f.readline().split()
            if int(fields[1]) != 4:
                continue
            ie = int(fields[0])
            ntags = int(fields[2])
            if ntags < 4:
                raise SystemExit("ERROR: MSH2 ELEMENT %d HAS NO PARTITION TAG" % ie)
            if ie < 1 or ie > nelem:
                raise SystemExit("ERROR: GMSH DID NOT PRESERVE ELEMENT NUMBERING")
            rank_of[ie] = int(fields[6]) - 1

    check_partition(rank_of, nelem, parts)
    return rank_of


def read_partition_file(path, nelem, parts):
    rank_of = [-1]
    with open(path, "r") as f:
        for line in f:
            if line.strip():
                rank_of.append(int(line.split()[0]))
    if len(rank_of) != nelem + 1:
        raise SystemExit("ERROR: PARTITION FILE HAS %d ENTRIES, CASE HAS %d TETS"
                         % (len(rank_of) - 1, nelem))
    check_partition(rank_of, nelem, parts)
    return rank_of


def check_partition(rank_of, nelem, parts):
    used = set()
    for ie in range(1, nelem + 1):
        r = rank_of[ie]
        if r < 0 or r >= parts:
            raise SystemExit("ERROR: ELEMENT %d HAS RANK %d OUTSIDE 0..%d" % (ie, r, parts - 1))
        used.add(r)
    if len(used) != parts:
        raise SystemExit("ERROR: ONLY %d OF %d RANKS RECEIVED ELEMENTS" % (len(used), parts))


# ================================
# BUILD AND WRITE THE RANK-LOCAL CASES
# ================================
def export_ranks(case, rank_of, parts, name, out_dir, par, matprop):
    nelem = case["nelem"]
    npoin = case["npoin"]
    elements = case["elements"]

    # ================================
    # STEP 1: NODE OWNER = LOWEST RANK TOUCHING THE NODE
    #         node_ranks[ip] = SET OF RANKS THAT HAVE THE NODE
    # ================================
    owner = [parts] * (npoin + 1)
    node_ranks = [None] * (npoin + 1)
    rank_elements = [[] for r in range(parts)]
    for ie in range(1, nelem + 1):
        r = rank_of[ie]
        rank_elements[r].append(ie)
        for ip in elements[ie]:
            if r < owner[ip]:
                owner[ip] = r
            if node_ranks[ip] is None:
                node_ranks[ip] = set()
            node_ranks[ip].add(r)
    for ip in range(1, npoin + 1):
        if node_ranks[ip] is None:
            raise SystemExit("ERROR: NODE %d IS NOT USED BY ANY TETRAHEDRON" % ip)

    # ================================
    # STEP 2: BOUNDARY FACES PER RANK, IN SERIAL ORDER
    # ================================
    rank_faces = [[] for r in range(parts)]
    for face in case["faces"]:
        rank_faces[rank_of[face[3]]].append(face)

    manifest_ranks = []
    for r in range(parts):
        # ================================
        # STEP 3: LOCAL NODE NUMBERING (OWNED FIRST, THEN GHOSTS)
        # ================================
        used = set()
        for ie in rank_elements[r]:
            used.update(elements[ie])
        owned = sorted(ip for ip in used if owner[ip] == r)
        ghost = sorted(ip for ip in used if owner[ip] != r)
        local_nodes = owned + ghost
        local_of = {}
        for il, ip in enumerate(local_nodes):
            local_of[ip] = il + 1

        local_elem_of = {}
        for il, ie in enumerate(rank_elements[r]):
            local_elem_of[ie] = il + 1

        # ================================
        # STEP 4: HALO LISTS
        # ================================
        send = {}
        recv = {}
        for ip in owned:
            for q in node_ranks[ip]:
                if q != r:
                    send.setdefault(q, []).append(local_of[ip])
        for ip in ghost:
            recv.setdefault(owner[ip], []).append(local_of[ip])
        neighbours = sorted(set(send) | set(recv))

        # ================================
        # STEP 5: WRITE THE FILES
        # ================================
        rank_id = "%04d" % r
        rank_dir = os.path.join(out_dir, "rank_" + rank_id)
        os.makedirs(rank_dir)
        base = os.path.join(rank_dir, name + "_rank_" + rank_id)
        faces = rank_faces[r]

        with open(base + ".plt", "w") as f:
            f.write("%10d%10d%10d\n" % (len(rank_elements[r]), len(local_nodes), len(faces)))
            for il, ie in enumerate(rank_elements[r]):
                n = [local_of[ip] for ip in elements[ie]]
                f.write("%10d%10d%10d%10d%10d\n" % (il + 1, n[0], n[1], n[2], n[3]))
            for il, ip in enumerate(local_nodes):
                f.write("%10d%s\n" % (il + 1, case["coord_text"][ip]))
            for face in faces:
                f.write("%10d%10d%10d%10d%10d\n" % (
                    local_of[face[0]], local_of[face[1]], local_of[face[2]],
                    local_elem_of[face[3]], face[4]))

        bc_counts = {}
        for face in faces:
            bc_counts[face[4]] = bc_counts.get(face[4], 0) + 1
        with open(base + ".bco", "w") as f:
            f.write("%10d%10d\n" % (len(bc_counts), case["mass_flow_flag"]))
            for bc in sorted(bc_counts):
                f.write("%10d%10d\n" % (bc, bc))

        material_counts = {}
        with open(base + ".material", "w") as f:
            for il, ie in enumerate(rank_elements[r]):
                n = [local_of[ip] for ip in elements[ie]]
                m = case["material"][ie]
                material_counts[m] = material_counts.get(m, 0) + 1
                f.write("%10d%10d%10d%10d%10d%6d\n" % (il + 1, n[0], n[1], n[2], n[3], m))

        shutil.copyfile(par, base + ".par")
        shutil.copyfile(matprop, base + ".matprop")

        with open(base + ".mpi", "w") as f:
            f.write("CBS3D_MPI_PARTITION 1\n")
            f.write("PARTITION %d %d %d\n" % (r + 1, r, parts))
            f.write("GLOBAL %d %d %d\n" % (nelem, npoin, case["nboun"]))
            f.write("LOCAL %d %d %d %d %d\n" % (len(rank_elements[r]), len(local_nodes),
                                                len(owned), len(ghost), len(faces)))
            f.write("NODES\n")
            for il, ip in enumerate(local_nodes):
                f.write("%d %d %d\n" % (il + 1, ip, owner[ip]))
            f.write("ELEMENTS\n")
            for il, ie in enumerate(rank_elements[r]):
                f.write("%d %d\n" % (il + 1, ie))
            f.write("NEIGHBOURS %d\n" % len(neighbours))
            for q in neighbours:
                s = send.get(q, [])
                v = recv.get(q, [])
                f.write("NEIGHBOUR %d %d %d\n" % (q, len(s), len(v)))
                f.write(" ".join(["SEND"] + [str(i) for i in s]) + "\n")
                f.write(" ".join(["RECV"] + [str(i) for i in v]) + "\n")
            f.write("END\n")

        manifest_ranks.append({
            "mpi_rank": r,
            "partition_id": r + 1,
            "case_base": "rank_%s/%s_rank_%s" % (rank_id, name, rank_id),
            "owned_tetrahedra": len(rank_elements[r]),
            "local_nodes": len(local_nodes),
            "owned_nodes": len(owned),
            "ghost_nodes": len(ghost),
            "physical_boundary_triangles": len(faces),
            "neighbour_ranks": neighbours,
            "bc_counts": dict((str(k), v) for k, v in sorted(bc_counts.items())),
            "material_counts": dict((str(k), v) for k, v in sorted(material_counts.items())),
        })
        print("  rank %4d : %8d tets %8d nodes (%d ghost) %7d faces, neighbours %s"
              % (r, len(rank_elements[r]), len(local_nodes), len(ghost), len(faces), neighbours))

    return manifest_ranks


# ================================
# MANIFEST AND COUNTS
# ================================
def write_manifest(out_dir, name, parts, case, ranks, par, matprop):
    owned = [x["owned_tetrahedra"] for x in ranks]
    owned_nodes = sum(x["owned_nodes"] for x in ranks)
    mean = float(sum(owned)) / parts
    manifest = {
        "format": "CBS3D rank-local text partition v1",
        "case_name": name,
        "parts": parts,
        "rank_mapping": "MPI rank = Gmsh partition ID - 1",
        "global": {"tetrahedra": case["nelem"], "nodes": case["npoin"],
                   "physical_boundary_triangles": case["nboun"]},
        "parameter_file": os.path.abspath(par),
        "material_property_file": os.path.abspath(matprop),
        "checks": {
            "owned_tetrahedra_sum_matches_global": sum(owned) == case["nelem"],
            "owned_nodes_sum_matches_global": owned_nodes == case["npoin"],
            "boundary_faces_sum_matches_global":
                sum(x["physical_boundary_triangles"] for x in ranks) == case["nboun"],
            "all_ranks_have_owned_tetrahedra": min(owned) > 0,
        },
        "load_balance": {
            "owned_tetrahedra_min": min(owned),
            "owned_tetrahedra_max": max(owned),
            "owned_tetrahedra_mean": mean,
            "owned_tetrahedra_max_over_mean": max(owned) / mean,
            "ghost_nodes_total": sum(x["ghost_nodes"] for x in ranks),
        },
        "ranks": ranks,
    }
    with open(os.path.join(out_dir, "partition_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")

    with open(os.path.join(out_dir, "partition_counts.csv"), "w") as f:
        f.write("partition_id,mpi_rank,owned_tetrahedra,local_nodes,owned_nodes,"
                "ghost_nodes,physical_boundary_triangles,neighbour_count,case_base\n")
        for x in ranks:
            f.write("%d,%d,%d,%d,%d,%d,%d,%d,%s\n" % (
                x["partition_id"], x["mpi_rank"], x["owned_tetrahedra"], x["local_nodes"],
                x["owned_nodes"], x["ghost_nodes"], x["physical_boundary_triangles"],
                len(x["neighbour_ranks"]), x["case_base"]))
    return manifest


# ================================
# MAIN
# ================================
def main():
    parser = argparse.ArgumentParser(
        description="Split a serial CBS3D case into rank-local MPI cases.")
    parser.add_argument("base", help="Serial case base path without extension")
    parser.add_argument("--parts", type=int, required=True, help="Number of MPI ranks")
    parser.add_argument("--out-dir", required=True, help="Partition root passed to cbs3d_parallel")
    parser.add_argument("--case-name", default=None, help="Rank case name (default: serial base name)")
    parser.add_argument("--par", default=None, help=".par for every rank (default: <base>.par)")
    parser.add_argument("--matprop", default=None, help=".matprop for every rank (default: <base>.matprop)")
    parser.add_argument("--partition-file", default=None,
                        help="One 0-based rank per tetrahedron instead of running Gmsh")
    parser.add_argument("--partition-msh", default=None,
                        help="Old-style partitioned MSH2 whose element tags are the CBS ids")
    parser.add_argument("--gmsh", default="gmsh", help="Gmsh executable (default: gmsh)")
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing --out-dir")
    args = parser.parse_args()

    base = args.base
    for ext in ("plt", "bco", "material"):
        if not os.path.isfile(base + "." + ext):
            raise SystemExit("ERROR: MISSING " + base + "." + ext)
    par = args.par or base + ".par"
    matprop = args.matprop or base + ".matprop"
    for path in (par, matprop):
        if not os.path.isfile(path):
            raise SystemExit("ERROR: MISSING " + path)
    if args.parts < 2:
        raise SystemExit("ERROR: --parts MUST BE AT LEAST 2 (USE THE SERIAL CASE FOR ONE RANK)")
    name = args.case_name or os.path.basename(base)

    if os.path.exists(args.out_dir):
        if not args.overwrite:
            raise SystemExit("ERROR: %s EXISTS; USE --overwrite TO REPLACE IT" % args.out_dir)
        shutil.rmtree(args.out_dir)
    os.makedirs(args.out_dir)

    print("READING SERIAL CASE " + base)
    case = read_serial_case(base)
    print("  %d tets, %d nodes, %d boundary faces" % (case["nelem"], case["npoin"], case["nboun"]))

    if args.partition_file:
        rank_of = read_partition_file(args.partition_file, case["nelem"], args.parts)
    elif args.partition_msh:
        rank_of = read_msh2_partition(args.partition_msh, case["nelem"], args.parts)
    else:
        work_dir = tempfile.mkdtemp(prefix="cbs_partition_")
        rank_of = partition_with_gmsh(case, args.parts, args.gmsh, work_dir)
        shutil.copyfile(os.path.join(work_dir, "gmsh_partition.log"),
                        os.path.join(args.out_dir, "gmsh_partition.log"))
        shutil.rmtree(work_dir)

    with open(os.path.join(args.out_dir, "element_partition.txt"), "w") as f:
        for ie in range(1, case["nelem"] + 1):
            f.write("%d\n" % rank_of[ie])

    print("WRITING %d RANKS TO %s" % (args.parts, args.out_dir))
    ranks = export_ranks(case, rank_of, args.parts, name, args.out_dir, par, matprop)
    manifest = write_manifest(args.out_dir, name, args.parts, case, ranks, par, matprop)

    lb = manifest["load_balance"]
    print("LOAD BALANCE max/mean = %.4f, ghost nodes = %d"
          % (lb["owned_tetrahedra_max_over_mean"], lb["ghost_nodes_total"]))
    failed = [k for k, v in manifest["checks"].items() if not v]
    if failed:
        print("RESULT: FAIL " + ", ".join(failed))
        return 1
    print("RESULT: PASS")
    print("RUN WITH: mpirun -np %d cbs3d_parallel %s %s" % (args.parts, name, args.out_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
