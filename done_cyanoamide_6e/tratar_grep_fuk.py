import os
import pandas as pd, numpy as np

def grepFukui():
    dataframes = []
    for fu in os.popen("ls 3_stage_rank_*_fk*.log").read().split("\n")[:3]:
        orbitais = os.popen("awk '/Summary of Natural Population Analysis:/{f=1;next} /=======================================================================/{f=0} f' " + fu).read()
        if orbitais.count("Natural Population") > 1:
            orbitais = np.array_split(orbitais.split("\n"), orbitais.count("Natural Population"))[-1][2:-1].tolist()
        else:
            orbitais = orbitais.split("\n")[3:-1]
        orbitais.pop(1)
        cabecalho = orbitais.pop(0).split()
        corpo = []
        for linha in orbitais:
            corpo.append(linha.split())
        print(cabecalho)
        dataframes.append(pd.DataFrame(corpo, columns=cabecalho))
    return dataframes
grepFukui()
#print(grepFukui()[1])