import meshio
import numpy
def temRepeticao(l):
    for x in l:
        if l.count(x) > 1:
            return True
    return False;


mesh = meshio.read('FisioPacer.vtk', 'vtk-ascii')
newTt = []

for cell in mesh.cells['tetra']:
    if not temRepeticao(list(cell)):
        newTt.append(cell)
print(len(newTt), " - ", len( mesh.cells['tetra']))
mesh.cells['tetra'] = numpy.array(newTt)
meshio.write_points_cells(
    "corrigido.vtk",
    mesh.points,
    mesh.cells,
    file_format='vtk-ascii',
    # Optionally provide extra data on points, cells, etc.
#    point_data=mesh.point_data,
 #   cell_data=mesh.cell_data,
  #  field_data=mesh.field_data
)

