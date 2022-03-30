%NProcShared=8
%Chk=opt1.chk
# AM1 Opt

 f5848555

0 1
C           0.92568        -0.01893         0.06417
C           2.44568        -0.01893         0.06417
N           3.70665        -0.01893         0.06417
H           0.55048        -0.80041         0.72806
H           0.55048        -0.20313        -0.94455
H           0.55048         0.94675         0.40901


--Link1--
%NProcShared=8
%Oldchk=opt1.chk
%Chk=opt2.chk
# m062x/6-311G(d,p) Opt Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 f5848555

0 1

--Link1--
%NProcShared=8
%Oldchk=opt2.chk
%Chk=fk+.chk
# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 f5848555

1 2

--Link1--
%NProcShared=8
%Oldchk=opt2.chk%Chk=fk-.chk
# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 f5848555

-1 2

