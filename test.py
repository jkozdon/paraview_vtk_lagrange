from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

pdi = self.GetInput()
pdo = self.GetOutput()

nelem = pdi.GetNumberOfCells()
newcells = vtk.vtkCellArray()
for i in range(nelem):
    cell = pdi.GetCell(i)
    N = cell.GetOrder(0)
    Nq = N + 1

    # corners
    newcells.InsertNextCell(Nq * Nq * Nq)
    newcells.InsertCellPoint(cell.GetPointId(0 + 0 * Nq + 0 * Nq * Nq))
    newcells.InsertCellPoint(cell.GetPointId(N + 0 * Nq + 0 * Nq * Nq))
    newcells.InsertCellPoint(cell.GetPointId(N + N * Nq + 0 * Nq * Nq))
    newcells.InsertCellPoint(cell.GetPointId(0 + N * Nq + 0 * Nq * Nq))
    newcells.InsertCellPoint(cell.GetPointId(0 + 0 * Nq + N * Nq * Nq))
    newcells.InsertCellPoint(cell.GetPointId(N + 0 * Nq + N * Nq * Nq))
    newcells.InsertCellPoint(cell.GetPointId(N + N * Nq + N * Nq * Nq))
    newcells.InsertCellPoint(cell.GetPointId(0 + N * Nq + N * Nq * Nq))

    # edges
    for n in range(1, N):
      i, j, k = n, 0, 0
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k = N, n, 0
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k = n, N, 0
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k =  0, n, 0
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))

    for n in range(1, N):
      i, j, k = n, 0, N
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k = N, n, N
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k = n, N, N
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k =  0, n, N
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))

    for n in range(1, N):
      i, j, k = 0, 0, n
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k = N, 0, n
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k = 0, N, n
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for n in range(1, N):
      i, j, k =  N, N, n
      newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))

    # faces
    for m in range(1, N):
        for n in range(1, N):
            i, j, k =  0, n, m
            newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for m in range(1, N):
        for n in range(1, N):
            i, j, k =  N, n, m
            newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for m in range(1, N):
        for n in range(1, N):
            i, j, k =  n, 0, m
            newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for m in range(1, N):
        for n in range(1, N):
            i, j, k =  n, N, m
            newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for m in range(1, N):
        for n in range(1, N):
            i, j, k =  n, m, 0
            newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))
    for m in range(1, N):
        for n in range(1, N):
            i, j, k =  n, m, N
            newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))

    # volume
    for k in range(1, N):
        for j in range(1, N):
            for i in range(1, N):
                newcells.InsertCellPoint(cell.GetPointId(i + j * Nq + k * Nq * Nq))

pdo.SetCells( 72, newcells)

# LGL points and barycentric weight (as computed by chebfun)
if Nq == 2:
    r = array([                    -1,                         1])
    w = array([                    -1,                         1])
elif Nq == 3:
    r = array([-1.000000000000000e+00,                         0,     1.000000000000000e+00,])
    w = array([-5.000000000000000e-01,     1.000000000000000e+00,    -5.000000000000000e-01,])
elif Nq == 4:
    r = array([-1.000000000000000e+00,    -4.472135954999579e-01,     4.472135954999579e-01,     1.000000000000000e+00,])
    w = array([-4.472135954999579e-01,     1.000000000000000e+00,    -1.000000000000000e+00,     4.472135954999579e-01,])
elif Nq == 5:
    r = array([-1.000000000000000e+00,    -6.546536707079771e-01,    -6.077163357286271e-64,     6.546536707079771e-01,     1.000000000000000e+00,])
    w = array([-3.749999999999999e-01,     8.749999999999999e-01,    -1.000000000000000e+00,     8.749999999999999e-01,    -3.749999999999999e-01,])
elif Nq == 6:
    r = array([-1.000000000000000e+00,    -7.650553239294646e-01,    -2.852315164806451e-01,     2.852315164806451e-01,     7.650553239294646e-01,     1.000000000000000e+00,])
    w = array([-3.466277250554179e-01,     8.259000647047562e-01,    -1.000000000000000e+00,     1.000000000000000e+00,    -8.259000647047562e-01,     3.466277250554179e-01,])
elif Nq == 7:
    r = array([-1.000000000000000e+00,    -8.302238962785669e-01,    -4.688487934707142e-01,                         0,     4.688487934707142e-01,     8.302238962785669e-01,     1.000000000000000e+00,])
    w = array([-3.124999999999999e-01,     7.534651069828725e-01,    -9.409651069828725e-01,     1.000000000000000e+00,    -9.409651069828725e-01,     7.534651069828725e-01,    -3.124999999999999e-01,])
elif Nq == 8:
    r = array([-1.000000000000000e+00,    -8.717401485096066e-01,    -5.917001814331423e-01,    -2.092992179024789e-01,     2.092992179024789e-01,     5.917001814331423e-01,     8.717401485096066e-01,     1.000000000000000e+00,])
    w = array([-2.942596405885329e-01,     7.147371237159649e-01,    -9.094210895526520e-01,     1.000000000000000e+00,    -1.000000000000000e+00,     9.094210895526520e-01,    -7.147371237159649e-01,     2.942596405885329e-01,])
elif Nq == 9:
    r = array([-1.000000000000000e+00,    -8.997579954114602e-01,    -6.771862795107377e-01,    -3.631174638261782e-01,     6.077163357286271e-64,     3.631174638261782e-01,     6.771862795107377e-01,     8.997579954114602e-01,     1.000000000000000e+00,])
    w = array([-2.734374999999998e-01,     6.674246433806348e-01,    -8.596291251131184e-01,     9.656419817324832e-01,    -1.000000000000000e+00,     9.656419817324832e-01,    -8.596291251131184e-01,     6.674246433806348e-01,    -2.734374999999998e-01,])
elif Nq == 10:
    r = array([-1.000000000000000e+00,    -9.195339081664589e-01,    -7.387738651055050e-01,    -4.779249498104445e-01,    -1.652789576663870e-01,     1.652789576663870e-01,     4.779249498104445e-01,     7.387738651055050e-01,     9.195339081664589e-01,     1.000000000000000e+00,])
    w = array([-2.604724104587959e-01,     6.379590749419513e-01,    -8.286142980093549e-01,     9.442590400325345e-01,    -1.000000000000000e+00,     1.000000000000000e+00,    -9.442590400325345e-01,     8.286142980093549e-01,    -6.379590749419513e-01,     2.604724104587959e-01,])

newPoints = vtk.vtkPoints()
x = empty([Nq, Nq, Nq])
y = empty([Nq, Nq, Nq])
z = empty([Nq, Nq, Nq])
tmp_x = empty(Nq)
tmp_y = empty(Nq)
tmp_z = empty(Nq)
h = 2/Nq
for e in range(nelem):
    for k in range(Nq):
        for j in range(Nq):
            for i in range(Nq):
                n = i + j * Nq + k * Nq * Nq + e * Nq * Nq * Nq
                coord = pdi.GetPoint(i + j * Nq + k * Nq * Nq + e * Nq * Nq * Nq)
                x[i, j, k], y[i, j, k], z[i, j, k] = coord[:3]
    # interpolate in x
    for k in range(Nq):
        for j in range(Nq):
            for i in range(Nq):
                re = i * 2.0 / N - 1
                numx = 0.0
                numy = 0.0
                numz = 0.0
                den = 0.0
                for n in range(Nq):
                    if re == r[n]:
                        numx = x[n, j, k]
                        numy = y[n, j, k]
                        numz = z[n, j, k]
                        den = 1
                        break
                    val = w[n] / (re - r[n])
                    numx = numx + x[n, j, k] * val
                    numy = numy + y[n, j, k] * val
                    numz = numz + z[n, j, k] * val
                    den = den + val
                tmp_x[i] = numx / den
                tmp_y[i] = numy / den
                tmp_z[i] = numz / den
            for i in range(Nq):
                x[i, j, k] = tmp_x[i]
                y[i, j, k] = tmp_y[i]
                z[i, j, k] = tmp_z[i]
    # interpolate in y
    for k in range(Nq):
        for i in range(Nq):
            for j in range(Nq):
                re = j * 2.0 / N - 1
                numx = 0.0
                numy = 0.0
                numz = 0.0
                den = 0.0
                for n in range(Nq):
                    if re == r[n]:
                        numx = x[i, n, k]
                        numy = y[i, n, k]
                        numz = z[i, n, k]
                        den = 1
                        break
                    val = w[n] / (re - r[n])
                    numx = numx + x[i, n, k] * val
                    numy = numy + y[i, n, k] * val
                    numz = numz + z[i, n, k] * val
                    den = den + val
                tmp_x[j] = numx / den
                tmp_y[j] = numy / den
                tmp_z[j] = numz / den
            for j in range(Nq):
                x[i, j, k] = tmp_x[j]
                y[i, j, k] = tmp_y[j]
                z[i, j, k] = tmp_z[j]
    # interpolate in z
    for i in range(Nq):
        for j in range(Nq):
            for k in range(Nq):
                re = k * 2.0 / N - 1
                numx = 0.0
                numy = 0.0
                numz = 0.0
                den = 0.0
                for n in range(Nq):
                    if re == r[n]:
                        numx = x[i, j, n]
                        numy = y[i, j, n]
                        numz = z[i, j, n]
                        den = 1
                        break
                    val = w[n] / (re - r[n])
                    numx = numx + x[i, j, n] * val
                    numy = numy + y[i, j, n] * val
                    numz = numz + z[i, j, n] * val
                    den = den + val
                tmp_x[k] = numx / den
                tmp_y[k] = numy / den
                tmp_z[k] = numz / den
            for k in range(Nq):
                x[i, j, k] = tmp_x[k]
                y[i, j, k] = tmp_y[k]
                z[i, j, k] = tmp_z[k]

    for k in range(Nq):
        for j in range(Nq):
            for i in range(Nq):
                n = i + j * Nq + k * Nq * Nq + e * Nq * Nq * Nq
                newPoints.InsertPoint(n, x[i, j, k], y[i, j, k], z[i, j, k])

pdo.SetPoints(newPoints)

np_di = dsa.WrapDataObject(pdi)
np_do = dsa.WrapDataObject(pdo)
for fldname in np_di.PointData.keys():
    fld = copy(np_di.PointData[fldname])
    fld[:] += 1
    np_do.PointData.append(fld, fldname)
