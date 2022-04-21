#!/usr/bin/python3
import argparse
from multiprocessing import Process, cpu_count
import numpy as np
import pandas as pd
import os
from uuid import uuid4


### Stagio 0 -> gera os pdb's
### Stagio 1 -> gera os com's com PRODUTO - REAGENTE
### Stagio 2 -> calcula fukui
def generate_conf(smiles, job_name):
  global args
  print("\033[1;93mGerando PDB\n\033[0;0m")
  os.system(f"obabel -:'{smiles}' -o pdb -O inicio.pdb --gen3d --ff GAFF -h -minimize")
  print("\033[1;93mGerando Conformeros\n\033[0;0m")
  os.system(f"obabel inicio.pdb -O conformeros.pdb --conformer --nconf 10 --writeconformers")
  conformeros = np.array_split(open("conformeros.pdb", "r").read().split("\n")[:-1], args["conf_num"])
  confor_index = 0
  print("\033[1;93mEscrevendo os Conformeros\n\033[0;0m")
  while (len(conformeros) > 0):
    conformero = conformeros.pop()
    with open(f"{args['storage_path']}/{job_name}/pdb/0_stage_conformero_{confor_index}.pdb", "w") as pdb:
      pdb.write("\n".join(conformero.tolist()))
    com_in = os.popen(f"obabel {args['storage_path']}/pdb/{job_name}/0_stage_conformero_{confor_index}.pdb -o com").readlines()[2:]
    with open(f"{args['storage_path']}/{job_name}/1_stage_conformero_{confor_index}.com", "w") as com:
      com.write(f"%NProcShared={args['threads']}\n")
      com.write(f"%Chk=1_stage_conformero_{confor_index}.chk\n")
      com.write("# AM1 Opt\n")
      com.write(f"\n{job_name}\n")
      com.write("".join(com_in))
      com.write("\n")
      com.write("--Link1--")
      com.write(f"%NProcShared={args['threads']}\n")
      com.write(f"%Oldchk=1_stage_conformero_{confor_index}.chk\n")
      com.write(f"%Chk=1_stage_conformero_{confor_index}_solv.chk\n")
      com.write("# m062x/6-311G(d,p) scrf=(SMD,solvent=water) scf=maxcycle=1000 maxdisk=200Gb\n")
      com.write(f"\n{job_name}\n")
      com.write("0 1")
    confor_index += 1
  print("\033[1;93mVá na pasta jobs e execute o 'run_all.py' para submeter todos jobs no slurm\n\033[0;0m")

def rankeamento():
  print("\033[1;93mRankeando os Conformeros (PRODUCT - REAGENT)\n\033[0;0m")
  log_gaus = os.popen("ls 1_stage_conformero_*.log").read().split()
  energias = {}
  while (len(log_gaus) > 0):
    nome_arquivo = log_gaus.pop()
    energias[nome_arquivo] = []
    grep = os.popen(f"grep 'HF=' {nome_arquivo}").split()
    for linhas in grep:
      for elemento in linhas.split("\\"):
        if "HF=" in elemento:
          energias[nome_arquivo].append(elemento)
    energias[nome_arquivo] = (float(energias[nome_arquivo][1][3:]) - float(energias[nome_arquivo][0][3:])) / 627.5
  rank = 0
  ordem_melhores = sorted(energias, key=energias.get)
  print("\033[1;93mEscrevendo o Rankeamento -> 1_stage_rank.log\n\033[0;0m")
  with open("1_stage_rank.log", "w") as rank:
    for key in ordem_melhores:
      rank.write("RANK NAME INDEX ENERGY(Kcal/mol) STATUS_FREQ \n\n")
      rank.write("{rank} {nome} {index} {energia}\n".format(rank = rank, index = int(key.replace('.','_').split('_')[-2]), nome = key, energia = energias[nome_arquivo]))    
      rank += 1
  exit(0)

def verifica_freq(nome_arquivo):
  if os.path.isfile(nome_arquivo) == True:
    frequencias = os.popen(f"grep 'Frequencies' {nome_arquivo}").read().split("\n")
  else :
    print("\033[1;93mArquivo de log não encontrado!.\n\033[0;0m")
    exit(1)
  matriz_verdade = []
  print("\033[1;93mEscrevendo Frequenica -> frequencia.log.\n\033[0;0m")
  with open(f"frequencia.log", "a") as log:
    log.write("{aste} Frequencia do {nome_arquivo} {aste}\n\n".format(aste=13*"#", nome_arquivo = nome_arquivo))
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
    log.write("{aste}\n".format(aste=49*"#"))
    log.write("{aste} Frequencia Inicial é {menor_freq} {aste}\n".format(aste=9*"#", menor_freq = matriz_verdade[0]))
    log.write("{aste} Frequencia Final é {maior_freq} {aste}\n".format(aste=10*"#", maior_freq = matriz_verdade[-1]))
    log.write("{aste}\n".format(aste=49*"#"))
    if matriz_verdade[0] >= 0:
      return True
    else:
      return False

## Preciso achar o melhor rank caso geral e criar 
def calcula_freq():
  ranks = pd.read_csv('1_stage_rank.log', sep=' ')
  ranks = ranks.set_index("RANK")
  atualiza_status = ranks.query("STATUS_FREQ == 'calculate").index
  if len(atualiza_status) > 0:
    print("\033[1;93mAtualizando 'log' Frequencia -> frequencia.log\n\033[0;0m")
    status_index = atualiza_status[0]
    if (verifica_freq(ranks.loc[status_index, "NAME"])):
      ranks.loc[status_index, "STATUS_FREQ"] = "Pass"
      ranks.to_csv("1_stage_rank.log", sep=' ')
      os.environ["STATUS_FREQ"] = "TRUE"
      exit(0)
    else:
      ranks.loc[status_index, "STATUS_FREQ"] = "Failed"
  print("\033[1;93mCalculando Frequencia\n\033[0;0m")
  rank_usado = ranks.loc[ranks['STATUS_FREQ'].isnull()].sort_values(by='RANK').index[0]
  nome = ranks.loc[rank_usado, 'NAME']
  with open(f"2_stage_rank_{rank_usado}.com", "w") as com:
    com.write(f"%NProcShared={os.environ['SLURM_NTASKS']}\n")
    com.write(f"%Oldchk=1_stage_conformero_{ranks.loc[rank_usado, 'INDEX']}.chk\n")
    com.write(f"%Chk=2_stage_rank_{rank_usado}.chk\n")
    com.write("# m062x/6-311G(d,p) freq=noraman Opt Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
    com.write(f" {nome}\n\n")
    com.write("0 1\n")
  ranks.iloc[rank_usado, "STATUS_FREQ"] = 'calculate'
  ranks.to_csv("1_stage_rank.log", sep=' ')
  exit(0)

def calcula_fukui():
  ranks = pd.read_csv('1_stage_rank.log', sep=' ')
  ligante = ranks.query("STATUS_FREQ == 'Pass")
  if len(ligante) == 1:
    print("\033[1;93mCalculando FUKUI\n\033[0;0m")
    with open(f"3_stage_rank_{ligante['RANK']}.com", "w") as com:
      com.write(f"%NProcShared={os.environ['SLURM_NTASKS']}\n")
      com.write(f"%Oldchk=2_stage_rank_{ligante['RANK']}.chk\n")
      com.write("%Chk=fk+.chk\n")
      com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
      com.write(f" {ligante['NAME']}\n\n")
      com.write("1 2\n\n")
      com.write("--Link1--\n")
      com.write(f"%NProcShared={os.environ['SLURM_NTASKS']}\n")
      com.write("%Oldchk=opt2.chk\n")
      com.write("%Chk=fk-.chk\n")
      com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
      com.write(f" {ligante['NAME']}\n\n")
      com.write("-1 2\n\n")
      com.write("--Link1--\n")
      com.write(f"%NProcShared={os.environ['SLURM_NTASKS']}\n")
      com.write("%Oldchk=opt2.chk\n")
      com.write("%Chk=fk0.chk\n")
      com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
      com.write(f" {ligante['NAME']}\n\n")
      com.write("0 2\n\n")
  elif len(ligante) > 1:
    print("Error --- Mais de um ligante")
  else:
    exit(0)

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
    the_file.write("export STATUS_FREQ=FALSE\n\n")
    for confor_index in range(args["conf_num"]):
      the_file.write(f"g09 < 1_stage_conformero_{confor_index}.com > 1_stage_conformero_{confor_index}.log\n\n")
    the_file.write("calcula_fukui.py --ranking 1")
    the_file.write("while [[ $STATUS_FREQ = 'FALSE' ]]")
    the_file.write("do")
    the_file.write("calcula_fukui.py --calc-freq 1")
    the_file.write("done")
    the_file.write("calcula_fukui.py --calc-fukui 1")
    the_file.write("date\n")
    the_file.write("echo \"End job\"\n\n")
  os.system(f"echo '{job_name}' >> {args['storage_path']}/jobs_index.txt")

def sub_rotina(smiles):
  global args
  while (len(smiles) > 0):
    job_name = str(uuid4())[:8]
    if not os.path.exists("{storage_path}/{job_name}".format(storage_path = args["storage_path"], job_name = job_name)):
      os.makedirs("{storage_path}/{job_name}".format(storage_path = args["storage_path"], job_name = job_name))
    generate_conf(smiles.pop(), job_name)
    cria_lanza(job_name)

def create_run_all():
  with open("jobs/run_all.py", 'w') as run_all:
    run_all.write("#!/usr/bin/python3\n")
    run_all.write("import os\n")
    run_all.write("pwd=os.environ['PWD']\n")
    run_all.write("pastas=os.popen('ls').read().split('\\n')\n")
    run_all.write("for pasta in pastas:\n")
    run_all.write("  os.chdir(f'{pwd}/{pasta}')\n")
    run_all.write("  os.system('sbatch lanzaFukui.sh')\n")

def run_jobs():
  global args
  os.chdir(f"{os.environ['PWD']}/{args['storage_path']}")
  os.system("./run_all.py")

def main():
  global args
  try:
    smiles = open(args["smiles_file"], "r").readlines()
  except:
    print("\033[1;93mVerifique se o arquivo existe!!!\033[0;81m")
    exit(1)
  if not os.path.exists("{storage_path}".format(storage_path = args["storage_path"])):
      os.makedirs("{storage_path}".format(storage_path = args["storage_path"]))
  create_run_all()
  smiles = np.array_split(smiles, cpu_count())
  processos = []
  for _ in range(cpu_count()):
    p = Process(target=sub_rotina, args=(smiles.pop().tolist()))
    p.start()
    processos.append(p)
  for p in processos:
    p.join()

def creditos():
  print("\033[1;35m{aste}".format(aste=('*'*69)))
  print("\033[1;35m{aste} \033[1;93mIdealizador  - Jerônimo Lameira \033[1;35m{aste}".format(aste=('*'*18)))
  print("\033[1;35m{aste} \033[1;93mCodificação  - Alfredo Gabriel  \033[1;35m{aste}".format(aste=('*'*18)))
  print("\033[1;35m{aste} \033[1;93mFluxo Lógica - Carlos Gabriel   \033[1;35m{aste}".format(aste=('*'*18)))
  print("\033[1;35m{aste} \033[1;93mFluxo Lógica - Clauber Henrique \033[1;35m{aste}".format(aste=('*'*18)))
  print("\033[1;35m{aste}".format(aste=('*'*69)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Otimiza com AM1 gaussian 09\n \033[1;37mContato: \033[1;34malfredogdso@gmail.com\033[0;0m")
    parser.add_argument("--smiles-file", "-sf", help="Arquivo '.smi' com as smiles.", type=str)
    parser.add_argument("--conf-num", "-cn", help="Quantidade de conformeros a serem usado \033[1;93m(Padrão 10)\033[0;0m", type=int)
    parser.add_argument("--threads", "-th", help="Quantidades de Threads \033[1;93m(Padrão 8)\033[0;0m", type=int)
    parser.add_argument("--time-job", "-tj", help="Tempo de execucao maxima no formato SLURM \033[1;93m(Padrão 5-00:00:00)\033[0;0m", type=str)
    parser.add_argument("--queue-job", "-qj", help="Fila de execucao \033[1;93m(Padrão 'cpu')\033[0;0m", type=str)
    parser.add_argument("--storage-path", "-st", help="Pasta para salvar os jobs", type=str)
    parser.add_argument("--extract-freq", "-ef", help="Digite o caminho para o arquivo a fim de extrai a Frequencia do '.log' \033[1;93m(Cria extrai a Frequencia para a pasta atual.)\033[0;0m", type=str)
    parser.add_argument("--ranking", "-r", help="Realiza o Rankemanto com base no (Produto - Reagente)/627.5", type=int)
    parser.add_argument("--calc-freq", "-cfreq", help="Calcula a Frequencia", type=int)
    parser.add_argument("--calc-fukui", "-cfukui", help="Calcula o índice de Fukui", type=int)
    parser.add_argument("--run-job", "-rj", help="Submete todos os jobs na pasta 1->Sim", type=int)
    args = {k: v for k, v in vars(parser.parse_args()).items()}
    if not args["extract_freq"] is None:
      verifica_freq(args["extract_freq"])
    if args["storage_path"] is None:
      args["storage_path"] = "jobs"
    elif args["storage_path"][-1] == "/":
      args["storage_path"] = args["storage_path"][:-1]
    if args["threads"] is None:
      args["threads"] = 8
    if args["queue_job"] is None:
      args["queue_job"] = "cpu"
    if args["time_job"] is None:
      args["time_job"] = "5-00:00:00"
    if args["run_job"] == True:
      run_jobs()
    creditos()