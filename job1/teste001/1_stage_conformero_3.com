%NProcShared=8
%Chk=1_stage_conformero_3.chk
# AM1 Opt

 job1/teste001/0_stage_conformero_3.pdb

0  1
C           9.28000         0.57200         2.56300
S           9.85400        -1.09500         2.37500
O          10.41500        -1.23000         1.04700
O          10.63600        -1.43200         3.55200
C           8.36400        -2.06200         2.43100
C           8.01300        -2.85900         1.34500
C           6.83200        -3.61100         1.39700
C           7.53100        -1.98700         3.54600
C           6.35600        -2.74700         3.58600
C           5.98000        -3.57100         2.51600
C           4.76900        -4.41900         2.56200
C           3.46000        -3.86800         2.01900
C           4.58300        -5.67400         3.03400
C           3.14700        -6.12700         2.88600
C           2.45800        -4.99000         2.16600
C           1.44300        -5.26600         1.07700
C           0.99000        -4.69500         2.38900
C           5.53300        -6.59600         3.65900
C           5.70700        -6.59400         5.05200
C           6.59300        -7.49200         5.65300
C           7.32100        -8.37100         4.86000
F           8.16700        -9.23400         5.43800
C           7.18200        -8.36800         3.47800
C           6.28300        -7.48200         2.87500
H          10.15300         1.22900         2.58500
H           8.73100         0.66000         3.50200
H           8.64400         0.83200         1.71400
H           8.63600        -2.90100         0.45300
H           6.58600        -4.21900         0.52900
H           7.78000        -1.34400         4.38700
H           5.74800        -2.68700         4.48300
H           3.57400        -3.58000         0.96800
H           3.12500        -3.00100         2.59800
H           2.69400        -6.28600         3.87200
H           3.08400        -7.05500         2.30800
H           1.18400        -6.29600         0.85700
H           1.42100        -4.61600         0.21000
H           0.42600        -5.34000         3.05400
H           0.66300        -3.66100         2.40200
H           5.13800        -5.90900         5.67400
H           6.71400        -7.51600         6.73100
H           7.76700        -9.05800         2.87900
H           6.17300        -7.48800         1.79300


--Link1--%NProcShared=8
%Oldchk=1_stage_conformero_3.chk
%Chk=1_stage_conformero_3_solv.chk
# m062x/6-311G(d,p) scrf=(SMD,solvent=water) scf=maxcycle=1000 maxdisk=200Gb

AUTHOR    GENERATED BY OPEN BABEL 2.3.2
0 1