%NProcShared=4
%Chk=exmp.chk
#n AM1 Opt

 test

0 1
O         -1.11080        0.92544        0.00000
H         -2.00480        0.94709       -0.00000
H         -0.26830        0.91881       -0.00000

--Link1--
%NProcShared=4
%OldChk=exmp.chk
%Chk=exmp1.chk
#n m062x/6-311(d,p) Opt geom=allcheck

 test

0 1

