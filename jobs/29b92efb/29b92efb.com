%NProcShared=8
%Chk=opt1.chk
#n AM1 Opt

 29b92efb


0  1
C           1.00676        -0.06684        -0.07408
C           2.52676        -0.06684        -0.07408
N           3.78773        -0.06684        -0.07408
H           0.63156        -0.87682         0.55472
H           0.63156        -0.20640        -1.08994
H           0.63156         0.88271         0.31299



0 1

--Link1--
%NProcShared=8
%Oldchk=opt1.chk
%Chk=opt2.chk
#n m062x/6-311G(d,p) Opt Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 29b92efb

0 1

--Link2--
%NProcShared=8
%Oldchk=opt2.chk
%Chk=fk+.chk
#n m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 29b92efb

1 2

--Link3--
%NProcShared=8
%oldchk=opt2.chk%Chk=fk-.chk
#n m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 29b92efb

-1 2

