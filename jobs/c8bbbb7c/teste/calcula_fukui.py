import argparse
from multiprocessing import Process, cpu_count
import numpy as np
import os
from uuid import uuid4


def cria_lanza(args, job_name):
  with open(f"{args['storage_path']}/{job_name}/lanzaFukui.sh", "w") as the_file:
    the_file.write("#!/bin/bash\n")
    the_file.write(f"#SBATCH -J {job_name}\n")
    the_file.write(f"#SBATCH --time {args['time_job']}\n")
    the_file.write(f"#SBATCH -p {args['queue_job']}\n")
    the_file.write(f"#SBATCH --ntasks {args['threads']}\n")
    the_file.write("\n")
    the_file.write("date; pwd;\n\n")
    the_file.write("module load gaussian/09\n\n")
    the_file.write("echo \"End job\"\n\n")
    the_file.write("date\n\n")
    the_file.write("module load gaussian/09\n\n")
    the_file.write(f"g09 < {job_name}.com > {job_name}.log\n\n")
    the_file.write("date\n")
  os.system(f"echo '{job_name}' >> {args['storage_path']}/jobs_index.txt")


def verifica_freq(nome_arquivo):
  frequencias = os.popen(f"grep 'Frequencies' {nome_arquivo}").read().split("\n")
  matriz_verdade = []
  with open(f"frequencia.log", "a") as log:
    #cabecalho = f"{30*"#"} Frequencia do {nome_arquivo} {30*"#"}\n"
    log.write("{aste} Frequencia do {nome_arquivo} {aste}\n".format(aste=30*"#", nome_arquivo = nome_arquivo))
    for linha in frequencias:
      colunas = linha.split()
      for coluna in colunas:
        if coluna != "--":
          log.write(f"{coluna} ")
      try:
        matriz_verdade.append(int(file[index]))
      except ValueError:
        pass
    matriz_verdade.sort()
    if matriz_verdade[0] >= 0:
      log.write("{aste}\n".format(aste=30*"#"))
      log.write("{aste} Maior Frequencia é {maior_freq} {aste}\n".format(aste=30*"#", maior_freq = matriz_verdade[-1]))
      log.write("{aste}\n".format(aste=30*"#"))
      exit(0)
    else:
      exit(1)


def cria_com(smiles, args):
  uid = str(uuid4())[:8]
  if not os.path.exists(args["storage_path"]):
    os.mkdir(args["storage_path"])
  smiles = smiles.tolist()
  while (len(smiles) != 0):
    smile = smiles.pop()
    if not os.path.exists(f"{args['storage_path']}/{uid}"):
      os.mkdir(f"{args['storage_path']}/{uid}")
    cordenadas = os.popen(f"obabel -:'{smile.split()[0]}' -o com --ff GAFF --gen3d -h --minimize").readlines()[6:]
    os.system(f"cp {os.getcwd()}/{__file__} {args['storage_path']}/{uid}/{uid}")
    with open(f"{args['storage_path']}/{uid}/{uid}_stage_1.com", "w") as com:
      com.write(f"%NProcShared={args['threads']}\n")
      com.write("%Chk=opt1.chk\n")
      com.write("# AM1 Opt freq=noraman\n\n")
      com.write(f" {uid}\n\n")
      com.write("0 1\n")
      com.writelines(cordenadas)
      com.write("\n")
    with open(f"{args['storage_path']}/{uid}/{uid}_stage_2.com", "w") as com:
      com.write(f"%NProcShared={args['threads']}\n")
      com.write("%Oldchk=opt1.chk\n")
      com.write("%Chk=opt2.chk\n")
      com.write("# m062x/6-311G(d,p) freq=noraman Opt Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
      com.write(f" {uid}\n\n")
      com.write("0 1\n\n")
    with open(f"{args['storage_path']}/{uid}/{uid}_stage_3.com", "w") as com:
      com.write(f"%NProcShared={args['threads']}\n")
      com.write("%Oldchk=opt2.chk\n")
      com.write("%Chk=fk+.chk\n")
      com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
      com.write(f" {uid}\n\n")
      com.write("1 2\n\n")
      com.write("--Link1--\n")
      com.write(f"%NProcShared={args['threads']}\n")
      com.write("%Oldchk=opt2.chk\n")
      com.write("%Chk=fk-.chk\n")
      com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
      com.write(f" {uid}\n\n")
      com.write("-1 2\n\n")
      com.write("--Link1--\n")
      com.write(f"%NProcShared={args['threads']}\n")
      com.write("%Oldchk=opt2.chk\n")
      com.write("%Chk=fk0.chk\n")
      com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
      com.write(f" {uid}\n\n")
      com.write("0 2\n\n")
    cria_lanza(args, uid)

def main(args):
  try:
    smiles = open(args["smiles_file"], "r").readlines()
  except:
    print("\033[1;93mVerifique se o arquivo existe!!!\033[0;81m")
    exit(1)
  smiles = np.array_split(smiles, cpu_count())
  processos = []
  for _ in range(cpu_count()):
    p = Process(target=cria_com, args=(smiles.pop(), args))
    p.start()
    processos.append(p)
  for p in processos:
    p.join()

def run_jobs():
  pastas = open("jobs/jobs_index.txt", 'r').readlines()
  for pasta in pastas:
    arquivos_na_pasta = os.popen(f'ls {pasta}/').read().split('\n')
    for arquivo in arquivos_na_pasta:
      if ".com" in arquivo:
        os.system(f"sbatch {pasta}/{arquivo}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Otimiza com AM1 gaussian 09\n \033[1;37mContato: \033[1;34malfredogdso@gmail.com\033[0;0m")
    parser.add_argument("--smiles-file", "-sf", help="Arquivo '.smi' com as smiles.", type=str)
    parser.add_argument("--storage-path", "-st", help="Pasta para salvar os jobs", type=str)
    parser.add_argument("--threads", "-th", help="Quantidades de Threads", type=int)
    parser.add_argument("--time-job", "-tj", help="Tempo de execucao maxima no formato SLURM", type=str)
    parser.add_argument("--queue-job", "-qj", help="Fila de execucao", type=str)
    parser.add_argument("--run-job", "-rj", help="Submete todos os jobs na pasta 1->Sim \033[1;93m(Este arquivo deve estar na pasta criada ao executar uma vez)\033[0;0m", type=int)
    parser.add_argument("--extract-freq", "-ef", help="Extrai a Frequencia do '.log' 1->Sim \033[1;93m(Cria extrai a Frequencia para a pasta atual.)\033[0;0m", type=int)
    args = {k: v for k, v in vars(parser.parse_args()).items()}
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
    else:
      if args["smiles_file"] is None:
        exit(1)
      main(args)