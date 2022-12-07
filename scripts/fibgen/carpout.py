import numpy as NP
def temRepeticao(l):
    for x in l:
        if l.count(x) > 1:
            return True
    return False
#
#
def eliminateDup(cells):
    saida = []
    for linha in cells['tetra']:
        if (not temRepeticao(list(linha))):
            saida.append(linha)
    return saida
 #
 #
def writefib(fib, sheet, normal):
    file = open("vent.fib", "w")
    file.write(str(len(fib))+"\n")
    for i in range(len(fib)):
        f, s, n = fib[i], sheet[i], normal[i]
        #fs=  NP.dot(NP.array(f), NP.array(s))
        #fn = NP.dot(NP.array(f), NP.array(n))
        #sn  = NP.dot(NP.array(s), NP.array(n))
        #if abs(fs)>1e-3  or abs(sn)>1e-3 or abs(fn)>1e-3:
        #    print(i,fs," ", fn," ",sn)
        file.write("%f %f %f %f %f %f %f %f %f\n" % (f[0], f[1], f[2], s[0], s[1], s[2], n[0], n[1], n[2]))
    file.close()
#
#
def writelem(mesh):
    file = open("vent.elem", "w")
    file.write("%d\n" % (len(mesh.cells['tetra'])))
    for it in mesh.cells['tetra']:
        file.write("Tt %d %d %d %d 1\n" % (it[0], it[1], it[2], it[3]))
    file.close()
#
#
def writepts(mesh):
    file = open("vent.pts", "w")
    file.write("%d\n" % (len(mesh.points)))
    for it in mesh.points:
        file.write("%f %f %f\n" % (it[0], it[1], it[2]))
    file.close()
#
#
def writeboud(mesh):
    fbound = open("vent.bound", "w")
    i = 0
    for it in mesh.points:
        if it[2] == 5.0:
            fbound.write('%d 0 0\n' % (i))
            fbound.write('%d 1 0\n' % (i))
            fbound.write('%d 2 0\n' % (i))
        i += 1
    fbound.close()
#
#
def writepress(mesh):
    fpress = open("vent.press", "w")
    i=0
    for it in mesh.cells['triangle']:
        if mesh.cell_data['triangle']['gmsh:geometrical'][i] == 30:
            fpress.write("%d %d %d\n" % (it[0], it[1], it[2]))
        i += 1
    fpress.close()
#
#
def writecarp(mesh, fibCell, sheetCell, normCell):
    writelem(mesh)
    writepts(mesh)
    writeboud(mesh)
    writepress(mesh)
    writefib(fibCell, sheetCell, normCell)