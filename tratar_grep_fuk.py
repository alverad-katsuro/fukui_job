from ast import match_case
from distutils.log import error
import os
import pandas as pd, numpy as np

def verifica_fuk_cargas(lista, value):
  passou = []
  for df in lista:
    if df.sum().Charge.round() == value:
      passou.append(df)
  return passou

def grepFukui():
    dataframes = {}
    arquivos = os.popen("ls 3_stage_rank_*_fk*.log").read().split("\n")[:3]
    arquivos.sort()
    for fu in arquivos:
        orbitais = os.popen("awk '/Summary of Natural Population Analysis:/{f=1;next} /=======================================================================/{f=0;next} f' " + fu).read()
        if orbitais.count("Natural Population") > 1:
            data_orbi = []
            np_arr = np.array_split(orbitais.split("\n")[1:], orbitais.count("Natural Population"))
            for orbi in np_arr:
                orbi = orbi.tolist()[2:-1]
                orbi.pop(1)
                cabecalho = orbi.pop(0).split()
                corpo = []
                for linha in orbi:
                    corpo.append(linha.split())
                df = pd.DataFrame(corpo, columns=cabecalho)
                df = df.apply(pd.to_numeric, errors='ignore')
                data_orbi.append(df)
            if "fk+" in fu:
                dataframes["fk+"] = (verifica_fuk_cargas(data_orbi, 1))
            elif "fk-" in fu:
                dataframes["fk-"] = (verifica_fuk_cargas(data_orbi, -1))
            elif "fk0" in fu:
                dataframes["fk0"] = (verifica_fuk_cargas(data_orbi, 0))
        else:
            print(fu)
            orbitais = orbitais.split("\n")[3:]
            orbitais.pop(1)
            cabecalho = orbitais.pop(0).split()
            corpo = []
            for linha in orbitais:
                corpo.append(linha.split())
            print(corpo)
            if "fk+" in fu:
                dataframes["fk+"] = pd.DataFrame(corpo[:-3], columns=cabecalho)
            elif "fk-" in fu:
                dataframes["fk-"] = pd.DataFrame(corpo[:-3], columns=cabecalho)
            elif "fk0" in fu:
                dataframes["fk0"] = pd.DataFrame(corpo[:-3], columns=cabecalho)
    return dataframes
grepFukui()



#print(grepFukui()[1])