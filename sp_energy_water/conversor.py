from genericpath import isfile
import os

class Elementos:
    ELEMENTOS = {1:"H", 2:"He", 3:"Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si",
15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni",
29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba",
57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At",
86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf",
99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds", 111: "Rg",
112: "Uub", 114: "Uuq"}

    def convert_string(string) -> str:
        if (len(string) > 23):
            array = string.split()
            array[0] = Elementos.ELEMENTOS[int(array[0])]
            return (f"{array[0]:3}{array[1]:>11} {array[2]:>10} {array[3]:>10}\n")
        else:
            return string


caminhos = [os.path.join("./arquivos_migracao", nome) for nome in os.listdir("./arquivos_migracao")]
arquivos = [arq for arq in caminhos if os.path.isfile(arq)]
xyz = [arq for arq in arquivos if arq.lower().endswith(".xyz")]
if (not os.path.isdir("saida_arquivo")):
    os.makedirs("./saida_arquivo")
for arquivo in xyz:
    file = open(arquivo, "r")
    linha_cont = 0
    linha_total = 0
    for linha in file.readlines():
        if (linha_cont == 0):
            linha_total = (int(linha) + 1)
            linha_cont += 1
        elif (linha_cont == 1):
            array = arquivo.split("/")
            array[1] = "saida_arquivo"
            array[2] = linha
            path = (array[0] +"/"+ array[1] + "/" + array[2].strip() + ".gjf")
            if (os.path.isfile(path)):
                saida = open(path, "w")
            else: 
                saida = open(path, "a")
            saida.write("%nprocshared=4\n")
            saida.write("%mem=256Mb\n")
            saida.write(("%chk=" + array[2].strip() +".chk "+ "\n"))
            saida.write("# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb \n")
            saida.write("alkyne ps \n")
            saida.write("0 1\n")
            linha_cont += 1
        else:
            if (linha_cont == linha_total):
                saida.write(Elementos.convert_string(linha))
                linha_cont = 0
                saida.close()
            else:
                saida.write(Elementos.convert_string(linha))
                linha_cont += 1

