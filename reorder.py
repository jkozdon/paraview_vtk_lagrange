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
