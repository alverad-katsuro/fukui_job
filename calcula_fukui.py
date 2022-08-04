#!/usr/bin/python3
import argparse, os
from multiprocessing import Process, cpu_count
import numpy as np
import pandas as pd
from uuid import uuid4


def generate_mol(mol, job_name):
  global args
  print("\033[1;34mGerando PDB\n\033[0;38m", flush=True)
  os.environ["OMP_NUM_THREADS"] = "1"
  print(mol)
  print(f"obabel -:'{mol}' -o pdb -O {args['storage_path']}/{job_name}/pdb/mol.pdb --gen3d --ff GAFF -h --minimize")
  os.system(f"obabel -:'{mol}' -o pdb -O {args['storage_path']}/{job_name}/pdb/mol.pdb --gen3d --ff GAFF -h --minimize")
  os.environ.pop("OMP_NUM_THREADS", None)
  print("\033[1;34mEscrevendo fukui job\n\033[0;38m", flush=True)
  com_in_mol = os.popen(f"obabel {args['storage_path']}/{job_name}/pdb/mol.pdb -o com").readlines()[4:]
  with open(f"{args['storage_path']}/{job_name}/1_stage_mol.com", "w") as com:
    com.write(f"%NProcShared={args['threads']}\n")
    com.write("%Chk=mol_opt.chk\n")
    com.write("#n m062x/6-311G(d,p) Opt\n")
    com.write(f"\n {job_name}\n")
    com.write("".join(com_in_mol))
  with open(f"{args['storage_path']}/{job_name}/2_stage_fukui0.com", "w") as com:
    com.write(f"%NProcShared={args['threads']}\n")
    com.write(f"%Oldchk=mol_opt.chk\n")
    com.write("%Chk=fk0.chk\n")
    com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
    com.write(f" {job_name}_fk+\n\n")
    com.write("0 1\n\n")
  with open(f"{args['storage_path']}/{job_name}/2_stage_fukui-.com", "w") as com:
    com.write(f"%NProcShared={args['threads']}\n")
    com.write(f"%Oldchk=mol_opt.chk\n")
    com.write("%Chk=fk-.chk\n")
    com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
    com.write(f" {job_name}_fk-\n\n")
    com.write("1 2\n\n")
  with open(f"{args['storage_path']}/{job_name}/2_stage_fukui+.com", "w") as com:
    com.write(f"%NProcShared={args['threads']}\n")
    com.write(f"%Oldchk=mol_opt.chk\n")
    com.write("%Chk=fk+.chk\n")
    com.write("# m062x/6-311G(d,p) SP Pop=NBO geom=check scf=maxcycle=1000 maxdisk=100Gb\n\n")
    com.write(f" {job_name}_fk0\n\n")
    com.write("-1 2\n\n")
  print("\033[1;32mVá na pasta jobs e execute o 'run_all.py' para submeter todos jobs no slurm\n\033[0;38m", flush=True)

def calcula_fukui():
  print(f"\033[1;34mStart 2_stage_fukui0.com\033[0;38m", flush=True)
  if "Normal termination" in os.popen(f"tail -n 1 2_stage_fukui0.log").read():
    print(f"\033[1;33m1_stage_mol.com já foi calculado\033[0;38m", flush=True)
  else:
    os.system(f"g09 < 2_stage_fukui0.com > 2_stage_fukui0.log")
  print(f"\033[1;34mStart 2_stage_fukui+.com\033[0;38m", flush=True)
  if "Normal termination" in os.popen(f"tail -n 1 2_stage_fukui+.log").read():
    print(f"\033[1;33m2_stage_fukui+.com já foi calculado\033[0;38m", flush=True)
  else:
    os.system(f"g09 < 2_stage_fukui+.com > 2_stage_fukui+.log")
  print(f"\033[1;34mStart 2_stage_fukui-.com\033[0;38m", flush=True)
  if "Normal termination" in os.popen(f"tail -n 1 2_stage_fukui-.log").read():
    print(f"\033[1;33m2_stage_fukui- já foi calculado\033[0;38m", flush=True)
  else:
    os.system(f"g09 < 2_stage_fukui-.com > 2_stage_fukui-.log")
  print(f"\033[1;34mEscrevendo os CSV com o resultado do FUKUI\033[0;38m", flush=True)
  mulliken()
  nbo_charge()

  

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
    generate_mol(row.MOL, job_name)
    cria_lanza(job_name)
    os.system(f"echo '{job_name}' >> {args['storage_path']}/jobs_index.txt")
    os.system(f"cp {__file__} {args['storage_path']}/{job_name}/run_{__file__.split('/')[-1]}")

def run_job():
  print(f"\033[1;34mStart 1_stage_mol.com\033[0;38m", flush=True)
  if "Normal termination" in os.popen(f"tail -n 1 1_stage_mol.log").read():
    print(f"\033[1;33m1_stage_mol.com já foi calculado\033[0;38m", flush=True)
  else:
    os.system(f"g09 < 1_stage_mol.com > 1_stage_mol.log")  
  calcula_fukui()
  print("\n\033[1;34m------- End Job -------\033[0;38m", flush=True)
  creditos()

def verifica_fuk_cargas(lista, value):
  passou = []
  for df in lista:
    if df.sum().Charge.round() == value:
      passou.append(df)
  return passou[0]

def nbo_charge():
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
    parser.add_argument("--threads", "-th", help="Quantidades de Threads \033[1;93m(Padrão 8)\033[0;0m", type=int)
    parser.add_argument("--time-job", "-tj", help="Tempo de execucao maxima no formato SLURM \033[1;93m(Padrão 5-00:00:00)\033[0;0m", type=str)
    parser.add_argument("--queue-job", "-qj", help="Fila de execucao \033[1;93m(Padrão 'cpu')\033[0;0m", type=str)
    parser.add_argument("--storage-path", "-sp", help="Pasta para salvar os jobs", type=str)
    parser.add_argument("--sub-jobs", "-sj", help="Submete todos os jobs na pasta \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    parser.add_argument("--run-job", "-rj", help="Descreve algumas flags \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    parser.add_argument("--help-full", "-hf", help="Descreve algumas flags \033[1;93m(1 -> Sim) OBS:Função Interna\033[0;0m", type=int)
    args = {k: v for k, v in vars(parser.parse_args()).items()}
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
    if args["queue_job"] is None:
      args["queue_job"] = "cpu"
    if args["time_job"] is None:
      args["time_job"] = "5-00:00:00"
    if args["help_full"] == True:
      help_full()
    if args["run_job"] == True:
      run_job()
    if args["sub_jobs"] == True:
      submete_jobs()
    if not args["smiles_file"] is None:
      main()
    creditos()
