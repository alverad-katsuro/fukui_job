#!/usr/bin/python3
import argparse, os
from multiprocessing import Process, cpu_count
import numpy as np
import pandas as pd
from uuid import uuid4


### Stagio 0 -> gera os pdb's
### Stagio 1 -> gera os com's com PRODUTO - REAGENTE
### Stagio 2 -> calcula fukui
def generate_conf(product, reagent, job_name):
  global args
  print("\033[1;34mGerando PDB\n\033[0;38m", flush=True)
  os.environ["OMP_NUM_THREADS"] = "1"
  os.system(f"obabel -:'{reagent}' -o pdb -O {args['storage_path']}/{job_name}/pdb/reagent.pdb --gen3d --ff GAFF -h -minimize")
  os.system(f"obabel -:'{product}' -o pdb -O {args['storage_path']}/{job_name}/pdb/product.pdb --gen3d --ff GAFF -h -minimize")
  print("\033[1;34mGerando Conformeros do produto\n\033[0;38m", flush=True)
  os.system(f"obabel {args['storage_path']}/{job_name}/pdb/product1.pdb -O {args['storage_path']}/{job_name}/pdb/conformeros_prod.pdb --conformer --nconf {args['conf_num']} --writeconformers")
  os.environ.pop("OMP_NUM_THREADS", None)
  conf_gen = len(os.popen(f"grep END {args['storage_path']}/{job_name}/pdb/conformeros_prod.pdb").readlines())
  args["conf_num"] = conf_gen
  conformeros = np.array_split(open(f"{args['storage_path']}/{job_name}/pdb/conformeros_prod.pdb", "r").read().split("\n")[:-1], args["conf_num"])
  confor_index = 0
  print("\033[1;34mEscrevendo os Conformeros\n\033[0;38m", flush=True)
  com_in_reagent = os.popen(f"obabel {args['storage_path']}/{job_name}/pdb/reagent1.pdb -o com").readlines()[4:]
  with open(f"{args['storage_path']}/{job_name}/1_stage_reagent.com", "w") as com:
    com.write(f"%NProcShared={args['threads']}\n")
    com.write("#n m062x/6-311G(d,p) Opt\n")
    com.write(f"\n {job_name}\n")
    com.write("".join(com_in_reagent))
  while (len(conformeros) > 0):
    conformero = conformeros.pop()
    with open(f"{args['storage_path']}/{job_name}/pdb/0_stage_conformero_{confor_index}.pdb", "w") as pdb:
      pdb.write("\n".join(conformero.tolist()))
    com_in = os.popen(f"obabel {args['storage_path']}/{job_name}/pdb/0_stage_conformero_{confor_index}.pdb -o com").readlines()[4:]
    with open(f"{args['storage_path']}/{job_name}/1_stage_conformero_{confor_index}.com", "w") as com:
      com.write(f"%NProcShared={args['threads']}\n")
      com.write(f"%Chk=1_stage_conformero_{confor_index}.chk\n")
      com.write("#n Opt m062x/6-311G(d,p) scrf=(SMD,solvent=water) scf=maxcycle=1000\n")
      com.write(f"\n {job_name}\n")
      com.write("".join(com_in))
      #com.write("--Link1--\n")
      #com.write(f"%NProcShared={args['threads']}\n")
      #com.write(f"%Oldchk=1_stage_conformero_{confor_index}.chk\n")
      #com.write("#n Opt m062x/6-311G(d,p) geom=check scrf=(SMD,solvent=water) scf=maxcycle=1000\n")
      #com.write(f"\n {job_name}_solv\n")
      #com.write("\n0  1\n")
    confor_index += 1
  print("\033[1;32mVá na pasta jobs e execute o 'run_all.py' para submeter todos jobs no slurm\n\033[0;38m", flush=True)



### verificar as subpastas erro de logica
def rankeamento():
  if os.path.isfile("1_stage_rank.log") == True:
    print(f"\033[1;33mRANK já existente, pulando etapa\033[0;38m", flush=True)
    os.environ["STATUS_FREQ"] = "FALSE"
  else:
    print("\033[1;34mRankeando os Conformeros (PRODUCT - REAGENT)\n\033[0;38m", flush=True)
    log_gaus = os.popen("ls 1_stage_conformero_*.log").read().split()
    log_gaus_reg = os.popen("ls 1_stage_reagent.log").read().split()
    if len(log_gaus_reg) == 0:
      log_gaus_reg = os.popen("ls 1_stage/1_stage_reagent.log").read().split()
    if len(log_gaus) == 0:
      log_gaus = os.popen("ls 1_stage/1_stage_conformero_*.log").read().split()
    energias = {}
    reagent = []
    for elemento in os.popen(f"grep -A 1 'HF=' 1_stage_reagent.log").read().split('\\'):
        if "HF=" in elemento:
          reagent.append(elemento.replace("\n ", "").replace(" ", ""))
    try:
      reagent = float(reagent[0][3:])
    except IndexError:
      print(f"\033[1;33mErro no rankeamento do {nome_arquivo}, reagent sem HF= \033[0;38m", flush=True)
      exit(1)
    while (len(log_gaus) > 0):
      nome_arquivo = log_gaus.pop()
      energias[nome_arquivo] = []
      grep = os.popen(f"grep -A 1 'HF=' {nome_arquivo}").read().split('\\')
      for elemento in grep:
        if "HF=" in elemento:
          energias[nome_arquivo].append(elemento.replace("\n ", "").replace(" ", ""))
      try:
        energias[nome_arquivo] = (float(energias[nome_arquivo][0][3:]) - reagent) / 627.5
      except IndexError:
        print(f"\033[1;33mErro no rankeamento do {nome_arquivo}, um dos calculos não rodou (product HF=)\033[0;38m", flush=True)
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
        print(key)
        rank.write("{rank} {nome} {index} {energia}\n".format(rank = rank_index, index = int(key.replace('.','_').split('_')[-2]), nome = key.split('/')[1], energia = energias[key]))    
        rank_index += 1

def modulo_verifica_freq(nome_do_arquivo):
  global args
  frequencias = os.popen(f"grep 'Frequencies' {nome_do_arquivo}").read().split("\n")[:-1]
  if len(frequencias) == 0:
    ranks = pd.read_csv('1_stage_rank.log', sep=' ').set_index("RANK")
    status = ranks.query("STATUS_FREQ == 'calculate'")
    if (status.empty):
      print("\033[1;33mNão há 'STATUS_FREQ = calculate' no arquivo de RANK. (Alguma etapa foi pulada)\n\033[0;38m", flush=True)
      return "ERROR"
    else:
      ranks.loc[status.index[0], "STATUS_FREQ"] = "ERROR"
      ranks.to_csv("1_stage_rank.log", sep=' ')
      return "ERROR"
  else:
    matriz_verdade = []
    print("\033[1;34mEscrevendo Frequenica -> frequencia.log.\n\033[0;38m", flush=True)
    with open(f"frequencia.log", "a") as log:
      log.write("{aste} Frequencia do {nome_arquivo} {aste}\n\n".format(aste=13*"#", nome_arquivo = nome_do_arquivo))
      log.write
      for linha in frequencias:
        colunas = linha.split()
        for coluna in colunas:
          if coluna != "--":
            try:
              log.write(f"{float(coluna):5.4f} ")
              matriz_verdade.append(float(coluna))
            except ValueError:
              log.write(f"{coluna} ")
        log.write("\n")
      matriz_verdade.sort()
      log.write("\n{aste}\n".format(aste=49*"#"))
      log.write("{aste} Frequencia Inicial é {menor_freq} {aste_2}\n".format(aste=10*"#", aste_2=10*"#", menor_freq = matriz_verdade[0]))
      log.write("{aste} Frequencia Final é {maior_freq} {aste_2}\n".format(aste=10*"#", aste_2=10*"#", maior_freq = matriz_verdade[-1]))
      log.write("{aste}\n".format(aste=49*"#"))
      return matriz_verdade[0]

def verifica_freq(nome_do_arquivo):
  global args
  #if os.path.exists("2_stage"):
  if os.path.isfile(nome_do_arquivo) == True:
    return modulo_verifica_freq(nome_do_arquivo)
  elif os.path.isfile(f"2_stage/{nome_do_arquivo}") == True:
    return modulo_verifica_freq(f"2_stage/{nome_do_arquivo}")
  else:
    print("\033[1;33mArquivo de log não encontrado!.\n\033[0;38m", flush=True)
    exit(1)

## Preciso achar o melhor rank caso geral e criar 
def calcula_freq():
  ranks = pd.read_csv('1_stage_rank.log', sep=' ').set_index("RANK")
  status_calculate = ranks.query("STATUS_FREQ == 'calculate'").index
  if (not ranks.query("STATUS_FREQ == 'Pass'").empty):
    print("\033[1;33mPulando calculo de Frequencia, já existe um que atende os requisitos\n\033[0;38m", flush=True)
    os.environ["STATUS_FREQ"] = "TRUE"
  elif (not status_calculate.empty):
    print("\033[1;34mAtualizando 'log' Frequencia -> frequencia.log\n\033[0;38m", flush=True)
    rank_usado = status_calculate[0]
    return_freq = verifica_freq(f"2_stage_rank_{rank_usado}.log")
    if return_freq == 'ERROR':
      pass
    elif return_freq >= 0:
      ranks.loc[rank_usado, "STATUS_FREQ"] = "Pass"
      ranks.to_csv("1_stage_rank.log", sep=' ')
      os.environ["STATUS_FREQ"] = "TRUE"
    else:
      ranks.loc[rank_usado, "STATUS_FREQ"] = "Failed"
      ranks.to_csv("1_stage_rank.log", sep=' ')
      print("\033[1;34mTodos possuem frequencia negativa\033[1;30m", flush=True)
  else:
    print("\033[1;34mCalculando Frequencia\n\033[0;38m", flush=True)
    rank_usado = ranks.loc[ranks['STATUS_FREQ'].isnull()].sort_values(by='RANK')
    if len(rank_usado) == 0:
      os.environ["STATUS_FREQ"] = "NO_FREQ"
      print("\033[1;34mTodos possuem frequencia negativa\033[1;30m", flush=True)
    else:
      rank_usado = rank_usado.index[0]
      nome = ranks.loc[rank_usado, 'NAME']
      with open(f"2_stage_rank_{rank_usado}.com", "w") as com:
        com.write(f"%NProcShared={os.environ['threads']}\n")
        com.write(f"%Oldchk=1_stage/1_stage_conformero_{ranks.loc[rank_usado, 'INDEX']}.chk\n")
        com.write(f"%Chk=2_stage_rank_{rank_usado}.chk\n")
        com.write("# m062x/6-311G(d,p) freq=noraman Opt geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
        com.write(f" {nome}\n\n")
        com.write("0 1\n")
      ranks.loc[rank_usado, "STATUS_FREQ"] = 'calculate'
      os.system(f"g09 < 2_stage_rank_{rank_usado}.com > 2_stage_rank_{rank_usado}.log")
      ranks.to_csv("1_stage_rank.log", sep=' ')

def cria_lanza(job_name):
  global args
  with open(f"{args['storage_path']}/{job_name}/lanzaFukui.sh", "w") as the_file:
    the_file.write("#!/bin/bash\n")
    the_file.write(f"#SBATCH -J {job_name}\n")
    the_file.write(f"#SBATCH --time {args['time_job']}\n")
    the_file.write(f"#SBATCH -p {args['queue_job']}\n")
    the_file.write(f"#SBATCH --ntasks {args['threads']}\n")
    the_file.write("\n")
    the_file.write("date; pwd;\n\n")
    the_file.write("module load gaussian/09\n\n")
    for env in args.keys():
      the_file.write(f"export {env}={args[env]}\n")
    the_file.write("chmod +x run_calcula_fukui.py\n")
    the_file.write("./run_calcula_fukui.py --run-job 1\n\n")
    the_file.write("date;")
    the_file.write("echo \"End job\"\n\n")

def create_run_all():
  global args
  with open(f"{args['storage_path']}/run_all.py", 'w') as run_all:
    run_all.write("#!/usr/bin/python3\n")
    run_all.write("import os\n")
    run_all.write("pwd=os.environ['PWD']\n")
    run_all.write("pastas=open('jobs_index.txt', 'r').read().split('\\n')[:-1]\n")
    run_all.write("for pasta in pastas:\n")
    run_all.write("  try:\n")
    run_all.write("    os.chdir(f'{pwd}/{pasta}')\n")
    run_all.write("    os.system('sbatch lanzaFukui.sh')\n")
    run_all.write("  except NotADirectoryError:\n")
    run_all.write("    pass\n")
  os.system(f"chmod +x {args['storage_path']}/run_all.py")

def submete_jobs():
  global args
  os.chdir(f"{os.environ['PWD']}/{args['storage_path']}")
  os.system("./run_all.py")

# dataframe -> NAME PRODUCT REAGENT
def sub_rotina(dataframe):
  global args
  for index, row in dataframe.iterrows():
    if (not dataframe.NAME.empty):
      job_name = f"{row.NAME}_{str(uuid4())[:2]}"
    else:
      job_name = f"{str(uuid4())[:8]}"
    if not os.path.exists("{storage_path}/{job_name}/pdb".format(storage_path = args["storage_path"], job_name = job_name)):
      os.makedirs("{storage_path}/{job_name}/pdb".format(storage_path = args["storage_path"], job_name = job_name))
    generate_conf(row.PRODUCT, row.REAGENT, job_name)
    cria_lanza(job_name)
    os.system(f"echo '{job_name}' >> {args['storage_path']}/jobs_index.txt")
    os.system(f"cp {__file__} {args['storage_path']}/{job_name}/run_{__file__.split('/')[-1]}")

def run_job():
  print(f"\033[1;34mStart 1_stage_reagent.com\033[0;38m", flush=True)
  if os.path.exists("1_stage"):
    if "Normal termination" in os.popen(f"tail -n 1 1_stage/1_stage_reagent.log").read():
      print(f"\033[1;33m1_stage_reagent.com já foi calculado\033[0;38m", flush=True)
    else:
      os.system(f"g09 < 1_stage/1_stage_reagent.com > 1_stage/1_stage_reagent.log")  
  else:
    if "Normal termination" in os.popen(f"tail -n 1 1_stage_reagent.log").read():
      print(f"\033[1;33m1_stage_reagent.com já foi calculado\033[0;38m", flush=True)
    else:
      if (os.system(f"g09 < 1_stage_reagent.com > 1_stage_reagent.log")) != 0:
        if ("Error termination request processed by link 9999" in os.popen(f"tail -n 20 1_stage_reagent.log").read()):
          print("\033[1;33mErro de base, substituindo por m062 e tentando novamente\033[0;38m", flush=True)
          os.system(f"g09 < 1_stage_reagent.com > 1_stage_reagent.log")
  for confor_index in range(int(os.environ["conf_num"])):
    print(f"\033[1;34mStart 1_stage_conformero_{confor_index}.com\033[0;38m", flush=True)
    if os.path.exists("1_stage"):
      if "Normal termination" in os.popen(f"tail -n 1 1_stage/1_stage_conformero_{confor_index}.log").read():
        print(f"\033[1;33m1_stage_conformero_{confor_index}.com já foi calculado\033[0;38m", flush=True)
      else:
        os.system(f"g09 < 1_stage/1_stage_conformero_{confor_index}.com > 1_stage/1_stage_conformero_{confor_index}.log")  
    else:
      if "Normal termination" in os.popen(f"tail -n 1 1_stage_conformero_{confor_index}.log").read():
        print(f"\033[1;33m1_stage_conformero_{confor_index}.com já foi calculado\033[0;38m", flush=True)
      else:
        if (os.system(f"g09 < 1_stage_conformero_{confor_index}.com > 1_stage_conformero_{confor_index}.log")) != 0:
          if ("Error termination request processed by link 9999" in os.popen(f"tail -n 20 1_stage_conformero_{confor_index}.log").read()):
            print(f"\033[1;33mErro no conformero {confor_index}\033[0;38m", flush=True)
  print("\033[1;34mCalculando RANK\033[0;38m", flush=True)
  rankeamento()
  if not os.path.exists("1_stage"):
      os.makedirs("1_stage")
  os.system("mv 1_stage_conformero_* 1_stage")
  while (os.environ["STATUS_FREQ"] == "FALSE"):
    calcula_freq()
  if not os.path.exists("2_stage"):
      os.makedirs("2_stage")
  os.system("mv 2_stage_rank_* 2_stage")
  if os.environ["STATUS_FREQ"] == "TRUE":
    print("\033[1;34mAchei a frequencia\033[0;38m", flush=True)
    print("\033[1;34mEnd Job!\033[0;38m", flush=True)
  creditos()

def main():
  global args
  try:
    dataframe = pd.read_csv(args["smiles_file"], sep=" ")
  except:
    print("\033[1;33mVerifique se o arquivo existe!!!\033[0;38m", flush=True)
    exit(1)
  if not os.path.exists("{storage_path}".format(storage_path = args["storage_path"])):
      os.makedirs("{storage_path}".format(storage_path = args["storage_path"]))
  dataframe = np.array_split(dataframe, cpu_count())
  processos = []
  for _ in range(cpu_count()):
    dataframe_pop = dataframe.pop()
    if len(dataframe_pop) > 0:
      p = Process(target=sub_rotina, args=([dataframe_pop]))
      p.start()
      processos.append(p)
  for p in processos:
    p.join()
  create_run_all()
  exit(0)

def creditos():
  print("\033[1;35m{aste}".format(aste=('*'*69)), flush=True)
  print("\033[1;35m{aste} \033[1;93mIdealizador  - Jerônimo Lameira \033[1;35m{aste}".format(aste=('*'*18)), flush=True)
  print("\033[1;35m{aste} \033[1;93mCodificação  - Alfredo Gabriel  \033[1;35m{aste}".format(aste=('*'*18)), flush=True)
  print("\033[1;35m{aste} \033[1;93mFluxo Lógica - Carlos Gabriel   \033[1;35m{aste}".format(aste=('*'*18)), flush=True)
  print("\033[1;35m{aste} \033[1;93mFluxo Lógica - Clauber Henrique \033[1;35m{aste}".format(aste=('*'*18)), flush=True)
  print("\033[1;35m{aste}".format(aste=('*'*69)), flush=True)
  exit(0)

def help_full():
  print("\033[1;32mA ser desenvolvido", flush=True)
  print("Funções marcadas com OBS são chamadas pelo slurm", flush=True)
  print("Só as use se ja tiver completado as etapas anteriores\033[0;38m", flush=True)
  exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Otimiza com AM1 gaussian 09\n \033[1;37mContato: \033[1;34malfredogdso@gmail.com\033[0;0m")
    parser.add_argument("--smiles-file", "-sf", help="Arquivo '.smi' com as smiles.", type=str)
    parser.add_argument("--conf-num", "-cn", help="Quantidade de conformeros a serem usado \033[1;93m(Padrão 10)\033[0;0m", type=int)
    parser.add_argument("--threads", "-th", help="Quantidades de Threads \033[1;93m(Padrão 8)\033[0;0m", type=int)
    parser.add_argument("--time-job", "-tj", help="Tempo de execucao maxima no formato SLURM \033[1;93m(Padrão 5-00:00:00)\033[0;0m", type=str)
    parser.add_argument("--queue-job", "-qj", help="Fila de execucao \033[1;93m(Padrão 'cpu')\033[0;0m", type=str)
    parser.add_argument("--storage-path", "-sp", help="Pasta para salvar os jobs", type=str)
    parser.add_argument("--extract-freq", "-ef", help="Digite o caminho para o arquivo a fim de extrai a Frequencia do '.log' \033[1;93m(Cria extrai a Frequencia para a pasta atual.)\033[0;0m", type=str)
    parser.add_argument("--ranking", "-r", help="Realiza o Rankemanto com base no (Produto - Reagente)/627.5 \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    parser.add_argument("--calc-freq", "-cfreq", help="Calcula a Frequencia \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    parser.add_argument("--sub-jobs", "-sj", help="Submete todos os jobs na pasta \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    parser.add_argument("--run-job", "-rj", help="Descreve algumas flags \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    parser.add_argument("--help-full", "-hf", help="Descreve algumas flags \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    args = {k: v for k, v in vars(parser.parse_args()).items()}
    if not args["extract_freq"] is None:
      verifica_freq(args["extract_freq"])
    if args["storage_path"] is None:
      args["storage_path"] = "jobs"
    elif args["storage_path"][0] == "/":
      args["storage_path"] = f"{args['storage_path']}"
    elif args["storage_path"][-1] == "/":
      args["storage_path"] = f"{os.environ['PWD']}/{args['storage_path'][:-1]}"
    else:
      args["storage_path"] = f"{os.environ['PWD']}/{args['storage_path']}"
    if args["threads"] is None:
      args["threads"] = 8
    if args["conf_num"] is None:
      args["conf_num"] = 10
    if args["queue_job"] is None:
      args["queue_job"] = "cpu"
    if args["time_job"] is None:
      args["time_job"] = "5-00:00:00"
    if args["calc_freq"] == True:
      calcula_freq()
    if args["ranking"] == True:
      rankeamento()
    if args["help_full"] == True:
      help_full()
    if args["run_job"] == True:
      run_job()
    if args["sub_jobs"] == True:
      submete_jobs()
    if not args["smiles_file"] is None:
      main()
    creditos()
