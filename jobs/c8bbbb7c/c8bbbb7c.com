%NProcShared=8
%Chk=opt1.chk
# AM1 Opt

 c8bbbb7c

0 1
C           1.05258        -0.02333        -0.08180
C           2.57258        -0.02333        -0.08180
N           3.83354        -0.02333        -0.08180
H           0.67737         0.94788        -0.41074
H           0.67737        -0.22407         0.92376
H           0.67737        -0.79380        -0.75843


--Link1--
%NProcShared=8
%Oldchk=opt1.chk
%Chk=opt2.chk
# m062x/6-311G(d,p) Opt Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 c8bbbb7c

0 1

--Link1--
%NProcShared=8
%Oldchk=opt2.chk
%Chk=fk+.chk
# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 c8bbbb7c

1 2

--Link1--
%NProcShared=8
%Oldchk=opt2.chk
%Chk=fk-.chk
# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 c8bbbb7c

-1 2

--Link1--
%NProcShared=8
%Oldchk=opt2.chk
%Chk=fk0.chk
# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 c8bbbb7c

0 2

