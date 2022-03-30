%nprocshared=4
%mem=512Mb
chk=opt1 
#n AM1 Opt

 name

0 1

--Link1--
%nprocshared=4
%mem=512Mb
%Oldchk=opt1.chk
%chk=opt2.chk 
#n b3lyp/6-311G(d,p) Opt Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 neutra 

0 1

--Link2--
%nprocshared=4
%mem=512Mb
%Oldchk=opt2.chk
%chk=opt3.chk
#n b3lyp/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 positiva 

1 2

--Link3--
%nprocshared=4
%mem=512Mb
%Oldchk=opt3.chk
%chk=opt4.chk
#n b3lyp/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb

 negativa 

-1 2

