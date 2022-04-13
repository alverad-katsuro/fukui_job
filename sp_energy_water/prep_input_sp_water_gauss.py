#!/usr/bin/python3
import os

arquivos = os.popen("ls *.log").read().split("\n")[:-1]
print(arquivos)
os.mkdir("agua_gauss_sg")
for arquivo in arquivos:
    nome = arquivo.split(".")[0] 
    try:
        com_in = os.popen(f"babel -i log -I {arquivo} -o com").readlines()[1:]
        with open(f"agua_gauss_sg/{nome}.com", "w") as com:
            com.write("%NProcShared=8\n")
            com.write("# m062x/6-311G(d,p) scrf=(SMD,solvent=water) scf=maxcycle=1000 maxdisk=200Gb\n")
            com.writelines(com_in)
        
    except:
        print(f"\n\nErro no {nome}\n\n")