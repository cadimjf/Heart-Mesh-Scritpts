import carpout
import meshio as MSH
import math as MATH
import numpy as NP

def interpolacao(z):
    return (-0.0031*z*z-0.0126*z+0.9878) -1
#
#
def getnormal(pt):
    x, y, z = pt
    d = MATH.sqrt(x * x + y * y)
    # se esta no apice, zera o vetor em XY e deixa paralelo a Z
    if x == 0 and y == 0:
        n = [0, 0, -1]
    else:
        zDir = interpolacao(z)
        n = [x / d, y / d, zDir]
    nn = NP.linalg.norm(n)
    n = n / nn
    return n
#
#
def getnormalpts(mesh):
    aNorms = []
    for pt in mesh.points:
        n = getnormal(pt)
        aNorms.append(n)
    return aNorms

#
#
def getbarycenter(pt1, pt2, pt3, pt4):
    bary = [0, 0, 0]
    for i in range(3):
        bary[i] = (pt1[i]+pt2[i]+pt3[i]+pt4[i])/4.0
    return bary
#
#
def getDistance(p1, p2):
    dif = p1[:]-p2[:]
    return MATH.sqrt((dif[0]**2) + (dif[1]**2) + (dif[2]**2))
##
def prepareSurfaces(mesh):
    endo=[False for j in range (len(mesh.points))]#30
    epi = [False for j in range (len(mesh.points))] #40
    i=0
    for elem in mesh.cells['triangle']:
        if mesh.cell_data['triangle']['gmsh:geometrical'][i] == 30:
            endo[elem[0]] = True
            endo[elem[1]] = True
            endo[elem[2]] = True
        if mesh.cell_data['triangle']['gmsh:geometrical'][i] == 40:
            epi[elem[0]] = True
            epi[elem[1]] = True
            epi[elem[2]] = True
        i+= 1
    return endo, epi
##
#
def closestPntOnSurface( point, surf, mesh):
    i = 0
    smallI = 0
    smallDiff = MATH.inf
    for pt in mesh.points:
        #se o ponto i pertence a superficie surf
        if surf[i]:
            d = getDistance(point, pt)
            if d < smallDiff:
                smallDiff = d
                smallI = i
        i = i + 1
    return smallI, smallDiff
#
#
def rotacao(v, a, u):
    x, y, z = v[0], v[1], v[2]
    cosa = MATH.cos(a)
    sena = MATH.sin(a)
    R = NP.array(
            [[cosa+x*x*(1-cosa),    x*y*(1-cosa)-z*sena,    x*z*(1-cosa)+y*sena],
             [y*x*(1-cosa)+z*sena,  cosa+y*y*(1-cosa),      y*z*(1-cosa)-x*sena],
             [z*x*(1-cosa)-y*sena, z*y*(1-cosa)+x*sena,     cosa+z*z*(1-cosa)]])

    return R.dot(NP.array(u))

def getfibersheet(pt, alphaEpi, alphaEndo, normal, endoSurf, epiSurf):
    # distancia do endocardio
    iEndo, dEndo = closestPntOnSurface(pt, endoSurf, mesh)
    # distancia do epicardio
    iEpi, dEpi = closestPntOnSurface(pt, epiSurf, mesh)
    wt = dEndo + dEpi
    alpha = alphaEpi * (dEpi / wt) + alphaEndo * (dEndo / wt)
    vet = [-normal[1], normal[0], 0]
    iFib = rotacao(normal, MATH.radians(alpha), vet)
    iSheet = NP.cross(normal, iFib)
    return iFib, iSheet
#
#
def getFiberPts(axis, mesh, endoSurf, epiSurf, alphaEpi, alphaEndo):
    aFib = []
    aSheet = []
    i=0
    for pt in mesh.points:
       iFib, iSheet = getfibersheet(pt, alphaEpi, alphaEndo, axis[i], endoSurf, epiSurf)
       aFib.append(iFib)
       aSheet.append(iSheet)
       i = i + 1
    return aFib, aSheet
#
#
def getAxisCells(mesh, alphaEpi, alphaEndo, endoSurf, epiSurf):
    aNorms = []
    aFib = []
    aSheet = []
    for cell in mesh.cells['tetra']:
        i1, i2, i3, i4 = cell[:]
        bary = getbarycenter(mesh.points[i1], mesh.points[i2], mesh.points[i3], mesh.points[i4])
        n = getnormal(bary)
        aNorms.append(n)
        iFib, iSheet = getfibersheet(bary, alphaEpi, alphaEndo, n, endoSurf, epiSurf)
        aFib.append(iFib)
        aSheet.append(iSheet)
    return aFib, aSheet, aNorms
#
#
def fibgen(mesh):
    alphaEpi = 90  # 90 degress
    alphaEndo = -90  # -90 degress
    endoSurf, epiSurf = prepareSurfaces(mesh)
    normPts = getnormalpts(mesh)
    fibPts, sheetPts = getFiberPts(normPts, mesh, endoSurf, epiSurf, alphaEpi, alphaEndo)
    fibCell, sheetCell, normCell = getAxisCells(mesh, alphaEpi, alphaEndo, endoSurf, epiSurf)
    return fibCell, sheetCell, normCell, fibPts, sheetPts, normPts
#
#
#
mesh = MSH.read("supercoarse.msh")
fibCell, sheetCell, normCell, fibPts, sheetPts, normPts = fibgen(mesh)

point_data = dict(
    fiberP=NP.array(fibPts, dtype=NP.float32),
    sheetP=NP.array(sheetPts, dtype=NP.float32),
    normalP=NP.array(normPts, dtype=NP.float32))
cell_data = dict(
    tetra=dict(
        fiberC=NP.array(fibCell, dtype=NP.float32),
        sheetC=NP.array(sheetCell, dtype=NP.float32),
        normalC=NP.array(normCell, dtype=NP.float32)
    )
)

#print(len(mesh.cells['tetra']))
mesh.cells['tetra'] = NP.array(carpout.eliminateDup(mesh.cells))
#print(len(mesh.cells['tetra']))
carpout.writecarp(mesh, fibCell, sheetCell, normCell)
del mesh.cells['triangle']
MSH.write_points_cells(
    "output.vtu",
    mesh.points,
    mesh.cells,
    file_format='vtu-ascii',
    # Optionally provide extra data on points, cells, etc.
     point_data= point_data,
     cell_data=cell_data
     #field_data=mesh.field_data
    )

