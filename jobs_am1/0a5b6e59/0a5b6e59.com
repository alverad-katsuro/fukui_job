%NProcShared=8
%Chk=opt1.chk
#n AM1 Opt

 0a5b6e59

0  1
C           0.94361        -0.04822        -0.06704
C           2.38761        -0.04822        -0.06704
N           3.53862        -0.04822        -0.06704
H           0.56818         0.76226        -0.69513
H           0.56818         0.09048         0.94890
H           0.56818        -0.99740        -0.45490



0 1

--Link1--
%NProcShared=8
%Oldchk=opt1.chk
%Chk=opt2.chk
#n m062x/6-311G(d,p) Opt Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 0a5b6e59

0 1

--Link2--
%NProcShared=8
%Oldchk=opt2.chk
%Chk=fk+.chk
#n m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 0a5b6e59

1 2

--Link3--
%NProcShared=8
%oldchk=opt2.chk%Chk=fk-.chk
#n m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 0a5b6e59

-1 2

