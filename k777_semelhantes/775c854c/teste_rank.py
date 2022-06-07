import os
def rankeamento():
  if os.path.isfile("1_stage_rank.log") == True:
    print(f"\033[1;33mRANK já existente, pulando etapa\033[0;38m", flush=True)
    os.environ["STATUS_FREQ"] = "FALSE"
  else:
    print("\033[1;34mRankeando os Conformeros (PRODUCT - REAGENT)\n\033[0;38m", flush=True)
    log_gaus = os.popen("ls 1_stage_conformero_*.log").read().split()
    if len(log_gaus) == 0:
      log_gaus = os.popen("ls 1_stage/1_stage_conformero_*.log").read().split()
    energias = {}
    while (len(log_gaus) > 0):
      nome_arquivo = log_gaus.pop()
      energias[nome_arquivo] = []
      grep = os.popen(f"grep -A 1 'HF=' {nome_arquivo}").read().split('\\')
      for elemento in grep:
        if "HF=" in elemento:
          energias[nome_arquivo].append(elemento.replace("\n ", "").replace(" ", ""))
      try:
        print(energias)
        energias[nome_arquivo] = (float(energias[nome_arquivo][1][3:]) - float(energias[nome_arquivo][0][3:])) / 627.5
      except IndexError:
        print(f"\033[1;33mErro no rankeamento do {nome_arquivo}, um dos calculos não rodou(product ou reagent)\033[0;38m", flush=True)
        energias.pop(nome_arquivo)
    if len(energias.keys()) == 0:
      os.environ["STATUS_FREQ"] = "ERROR_RANK"
      exit(1)
    else:
      os.environ["STATUS_FREQ"] = "FALSE"
    rank_index = 0
    ordem_melhores = sorted(energias, key=energias.get)
    print("\033[1;34mEscrevendo o Rankeamento -> 1_stage_rank.log\n\033[0;38m", flush=True)
    with open("1_stage_rank.log", "w") as rank:
      rank.write("RANK NAME INDEX ENERGY(Kcal/mol) STATUS_FREQ STATUS_FUKUI\n")
      for key in ordem_melhores:
        rank.write("{rank} {nome} {index} {energia}\n".format(rank = rank_index, index = int(key.replace('.','_').split('_')[-2]), nome = key, energia = energias[nome_arquivo]))    
        rank_index += 1

rankeamento()