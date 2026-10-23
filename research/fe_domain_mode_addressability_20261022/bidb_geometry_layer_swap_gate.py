#!/usr/bin/env python3
"""Source-structure-only Seitz/site matching for an external ASE/SQLite BiDB.

Use: python bidb_geometry_layer_swap_gate.py --db /path/to/bidb.db \
       --ids 980 975 2684 2685 8055

Reads but never copies/uploads the DB. Assumes the printed primitive
hexagonal lattice convention a1=(a,0), a2=(-a/2,sqrt(3)*a/2).
Operations are C3z, chosen vertical mirror mv, and pure horizontal
mirror mz about the slab midplane. Atom matching is bijective, same Z,
modulo INTEGER in-plane Bravais translations. Numeric tolerance accounts
for relaxation/roundoff (default 0.025 A).

This tests *STRUCTURAL* symmetry only. It does not validate magnetic
space groups, SOC, tight-binding orbitals, or Kubo responses. A failed
pure mz does NOT prove absence of all layer-swap Seitz operations.
"""
import argparse
import json
import sqlite3
import numpy as np

C3F=np.array([[0,-1],[1,-1]],dtype=int)
MVF=np.array([[-1,1],[0,1]],dtype=int)


def validate(conn, db_id, op, tolerance=0.025):
    row=conn.execute(
      "select numbers,positions,cell from systems where id=?",
      (db_id,)).fetchone()
    if row is None: raise ValueError("unknown database id "+str(db_id))
    numbers, positions, cell=row
    species=np.frombuffer(numbers,dtype=np.int32)
    coords=np.frombuffer(positions,dtype=np.float64).reshape(-1,3)
    lattice=np.frombuffer(cell,dtype=np.float64).reshape(3,3)
    f=coords@np.linalg.inv(lattice)
    zmid=(np.min(coords[:,2])+np.max(coords[:,2]))/2
    used=set();perm=[];worst=0.0
    for i in range(len(species)):
        p=f[i].copy()
        if op=="C3z":p[:2]=C3F@p[:2]
        elif op=="mv":p[:2]=MVF@p[:2]
        elif op=="mz":p[2]=(2*zmid-coords[i,2])/lattice[2,2]
        else:raise ValueError(op)
        matches=[]
        for j in range(len(species)):
            if species[j]!=species[i] or j in used:continue
            delta=p-f[j]
            delta[:2]-=np.rint(delta[:2])
            matches.append((np.linalg.norm(delta@lattice),j))
        if not matches:
            return {"pass":False,"matched":len(perm),"reason":"no species match"}
        error,j=min(matches)
        if error>tolerance:
            return {"pass":False,"matched":len(perm),
                    "first_unmatched_error_angstrom":round(float(error),6)}
        used.add(j);perm.append(j);worst=max(worst,error)
    assert len(set(perm))==len(species)
    return {"pass":True,"max_error_angstrom":round(float(worst),6),
            "permutation_zero_based":perm}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--db",required=True)
    p.add_argument("--ids",type=int,nargs="+",required=True)
    p.add_argument("--tolerance",type=float,default=0.025)
    a=p.parse_args()
    conn=sqlite3.connect("file:"+a.db+"?mode=ro",uri=True)
    conn.execute("pragma query_only=ON")
    for id in a.ids:
        print(json.dumps({"db_id":id,
          "operators":{op:validate(conn,id,op,a.tolerance)
                       for op in ("C3z","mv","mz")}},ensure_ascii=False))
    conn.close()


if __name__=="__main__":
    main()
