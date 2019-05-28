def interpolate(Nq, Nqk, r, w, fld, tmp, nelem):
    fld = reshape(fld, (Nq, Nq, Nqk, nelem), order='F')
    h = 2/Nq
    for e in range(nelem):
        # interpolate in r
        for k in range(Nqk):
            for j in range(Nq):
                for i in range(Nq):
                    re = i * 2.0 / (Nq - 1) - 1
                    numx = 0.0
                    den = 0.0
                    for n in range(Nq):
                        if re == r[n]:
                            numx = fld[n, j, k, e]
                            den = 1
                            break
                        val = w[n] / (re - r[n])
                        numx = numx + fld[n, j, k, e] * val
                        den = den + val
                    tmp[i] = numx / den
                for i in range(Nq):
                    fld[i, j, k, e] = tmp[i]
        # interpolate in s
        for k in range(Nqk):
            for i in range(Nq):
                for j in range(Nq):
                    re = j * 2.0 / (Nq - 1) - 1
                    numx = 0.0
                    den = 0.0
                    for n in range(Nq):
                        if re == r[n]:
                            numx = fld[i, n, k, e]
                            den = 1
                            break
                        val = w[n] / (re - r[n])
                        numx = numx + fld[i, n, k, e] * val
                        den = den + val
                    tmp[j] = numx / den
                for j in range(Nq):
                    fld[i, j, k, e] = tmp[j]
        # interpolate in t
        if Nqk > 1:
            for i in range(Nq):
                for j in range(Nq):
                    for k in range(Nq):
                        re = k * 2.0 / (Nq - 1) - 1
                        numx = 0.0
                        den = 0.0
                        for n in range(Nq):
                            if re == r[n]:
                                numx = fld[i, j, n, e]
                                den = 1
                                break
                            val = w[n] / (re - r[n])
                            numx = numx + fld[i, j, n, e] * val
                            den = den + val
                        tmp[k] = numx / den
                    for k in range(Nq):
                        fld[i, j, k, e] = tmp[k]

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
x = empty([Nq, Nq, Nqk, 1])
y = empty([Nq, Nq, Nqk, 1])
z = empty([Nq, Nq, Nqk, 1])
tmp = empty(Nq)
print("interp xyz")
for e in range(nelem):
    for k in range(Nqk):
        for j in range(Nq):
            for i in range(Nq):
                n = i + j * Nq + k * Nqk * Nq + e * Nqk * Nq * Nq
                coord = pdi.GetPoint(i + j * Nq + k * Nqk * Nq + e * Nqk * Nq * Nq)
                x[i, j, k, 0], y[i, j, k, 0], z[i, j, k, 0] = coord[:3]
    interpolate(Nq, Nqk, r, w, x, tmp, 1)
    interpolate(Nq, Nqk, r, w, y, tmp, 1)
    interpolate(Nq, Nqk, r, w, z, tmp, 1)

    for k in range(Nqk):
        for j in range(Nq):
            for i in range(Nq):
                n = i + j * Nq + k * Nqk * Nq + e * Nqk * Nq * Nq
                newPoints.InsertPoint(n, x[i, j, k, 0], y[i, j, k, 0], z[i, j, k, 0])

pdo.SetPoints(newPoints)

np_di = dsa.WrapDataObject(pdi)
np_do = dsa.WrapDataObject(pdo)
for fldname in np_di.PointData.keys():
    print(fldname)
    fld = copy(np_di.PointData[fldname])
    interpolate(Nq, Nqk, r, w, fld, tmp, size(fld) / (Nq * Nq * Nqk))
    np_do.PointData.append(fld, fldname)
