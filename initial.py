from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

pdi = self.GetInput()
pdo = self.GetOutput()
nelem = pdi.GetNumberOfCells()
Nq = pdi.GetCell(0).GetOrder(0)+1
