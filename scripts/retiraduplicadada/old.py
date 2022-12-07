fin = open("vent_duplicados.elem", "r")
saida = []
contents = fin.readlines()
for linha in contents:
    listLinha = linha.split()[1:5]

    if (not temRepeticao(listLinha)):
        saida.append(listLinha)

ffib = open("vent.fib", "w")
ffib.write(str(len(saida)) + "\n")

fout = open("vent.elem", "w")
fout.write(str(len(saida)) + "\n")
for linha in saida:
    fout.write('Tt %s %s %s %s 1\n' % (linha[0], linha[1], linha[2], linha[3]))
    ffib.write('1 0 0 0 1 0 0 0 1\n')
fout.close()
ffib.close()

fpts = open("vent.pts", "r")
contpts = fpts.readlines()
boud = []
cont = 0
for linha in contpts:
    listLinha = linha.split()

    if (len(listLinha) < 3):
        continue
    # print(listLinha[2], listLinha[2]=='5')
    if str(listLinha[2]) == '5':
        # está na base, condição de contorno
        boud.append([cont, 0, 0])
        boud.append([cont, 1, 0])
        boud.append([cont, 2, 0])
    cont = cont + 1
fpts.close()
fbound = open("vent.bound", "w")
fbound.write(str(len(boud)) + "\n")
for linha in boud:
    fbound.write('%s %s %s\n' % (linha[0], linha[1], linha[2]))
fbound.close()

import meshio

press = meshio.read('press.vtu', 'vtu-ascii')
ids = press.point_data['Ids']
fpress = open("vent.press", "w")
for tri in press.cells['triangle']:
    fpress.write('%d %d %d\n' % (ids[tri[0]], ids[tri[1]], ids[tri[2]]))
fpress.close()
