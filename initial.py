from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

pdi = self.GetInput()
pdo = self.GetOutput()
nelem = pdi.GetNumberOfCells()
Nq = pdi.GetCell(0).GetOrder(0)+1
if pdi.GetCell(0).GetNumberOfPoints() == Nq * Nq:
    Nqk = 1
elif pdi.GetCell(0).GetNumberOfPoints() == Nq * Nq * Nq:
    Nqk = Nq
else:
    raise Exception('cannot determine dimension')
