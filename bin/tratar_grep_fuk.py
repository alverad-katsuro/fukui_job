import os
import pandas as pd, numpy as np

def verifica_fuk_cargas(lista, value):
  passou = []
  for df in lista:
    if df.sum().Charge.round() == value:
      passou.append(df)
  return passou[0]

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
                dataframes["n-1"] = (verifica_fuk_cargas(data_orbi, 1))
            elif "fk-" in fu:
                dataframes["n+1"] = (verifica_fuk_cargas(data_orbi, -1))
            elif "fk0" in fu:
                dataframes["n"] = (verifica_fuk_cargas(data_orbi, 0))
        else:
            orbitais = orbitais.split("\n")[3:]
            orbitais.pop(1)
            cabecalho = orbitais.pop(0).split()
            corpo = []
            for linha in orbitais:
                corpo.append(linha.split())
            if "fk+" in fu:
                dataframes["n-1"] = pd.DataFrame(corpo[:-3], columns=cabecalho)
            elif "fk-" in fu:
                dataframes["n+1"] = pd.DataFrame(corpo[:-3], columns=cabecalho)
            elif "fk0" in fu:
                dataframes["n"] = pd.DataFrame(corpo[:-3], columns=cabecalho)
    fukui = pd.DataFrame(columns=["No", "Atom", "n", "n+1", "n-1", "fk+", "fk-", "fk0"])
    fukui["No"] = dataframes["n"].No
    fukui["Atom"] = dataframes["n"].Atom
    fukui["n"] = dataframes["n"].Charge
    fukui["n+1"] = dataframes["n+1"].Charge
    fukui["n-1"] = dataframes["n-1"].Charge
    fukui = fukui.apply(pd.to_numeric, errors='ignore')
    fukui["fk+"] = fukui["n+1"] - fukui["n"]
    fukui["fk-"] = fukui["n"] - fukui["n-1"]
    fukui["fk0"] = (fukui["n+1"] - fukui["n-1"]) / 2
    fukui.to_csv("fukui_nbo_charge.csv", index=False)
    return fukui

def mulliken():
    dataframes = {}
    arquivos = os.popen("ls 3_stage_rank_*_fk*.log").read().split("\n")[:3]
    arquivos.sort()
    for fu in arquivos:
        if "fk0" in fu:
            orbitais = os.popen("awk '/ Mulliken charges:/{f=1;next} / Sum of Mulliken charges = /{f=0;next} f' " + fu).read()
            orbitais = orbitais.split("\n")[1:-1]
            corpo = []
            for linha in orbitais:
                corpo.append(linha.split())
        else:
            orbitais = os.popen("awk '/ Mulliken charges and spin densities:/{f=1;next} / Sum of Mulliken charges = /{f=0;next} f' " + fu).read()
            orbitais = orbitais.split("\n")[1:-1]
            corpo = []
            for linha in orbitais:
                corpo.append(linha.split()[:-1])
        if "fk+" in fu:
            dataframes["n-1"] = pd.DataFrame(data=corpo, columns=['No', 'Atom','Charge'])
        elif "fk-" in fu:
            dataframes["n+1"] = pd.DataFrame(data=corpo, columns=['No', 'Atom','Charge'])
        elif "fk0" in fu:
            dataframes["n"] = pd.DataFrame(data=corpo, columns=['No', 'Atom','Charge'])

    fukui = pd.DataFrame(columns=["No", "Atom","n", "n+1", "n-1", "fk+", "fk-", "fk0"])
    fukui["No"] = dataframes["n"].No
    fukui["Atom"] = dataframes["n"].Atom
    fukui["n"] = dataframes["n"].Charge
    fukui["n+1"] = dataframes["n+1"].Charge
    fukui["n-1"] = dataframes["n-1"].Charge
    fukui = fukui.apply(pd.to_numeric, errors='ignore')
    fukui["fk+"] = fukui["n+1"] - fukui["n"]
    fukui["fk-"] = fukui["n"] - fukui["n-1"]
    fukui["fk0"] = (fukui["n+1"] - fukui["n-1"]) / 2
    fukui.to_csv("fukui_mulliken.csv", index=False)
    return fukui
grepFukui()
mulliken()
