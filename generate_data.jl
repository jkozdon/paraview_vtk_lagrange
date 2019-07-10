using Printf
try
  GaussQuadrature.neither
catch
  include("quadrature.jl")
end

# From the VTK source:
# https://gitlab.kitware.com/nick.laurenson/vtk/blob/d8aa5b89c622cf04b3322112fda420afc9d2a16d/Common/DataModel/vtkLagrangeHexahedron.cxx#L729-795
function PointIndexFromIJK(i, j, k, N)
  ibdy = (i == 0 || i == N)
  jbdy = (j == 0 || j == N)
  kbdy = (k == 0 || k == N)
  # How many boundaries do we lie on at once?
  nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0)

  if (nbdy == 3) # Vertex DOF
    return (i != 0 ? (j != 0 ? 2 : 1) : (j != 0 ? 3 : 0)) + (k != 0 ? 4 : 0)
  end

  offset = 8
  if (nbdy == 2) # Edge DOF
    if (!ibdy)
      # On i axis
      return (i - 1) +
      (j != 0 ? N - 1 + N - 1 : 0) +
      (k != 0 ? 2 * (N - 1 + N - 1) : 0) +
      offset
    end
    if (!jbdy)
      # On j axis
      return (j - 1) +
      (i != 0 ? N - 1 : 2 * (N - 1) + N - 1) +
      (k != 0 ? 2 * (N - 1 + N - 1) : 0) +
      offset
    end
    # !kbdy, On k axis
    offset += 4 * (N - 1) + 4 * (N - 1)
    return (k - 1) + (N - 1) * (i != 0 ? (j != 0 ? 3 : 1) : (j != 0 ? 2 : 0)) + offset
  end

  offset += 4 * (N - 1 + N - 1 + N - 1)
  if (nbdy == 1) # Face DOF
    if (ibdy) # On i-normal face
      return (j - 1) + ((N - 1) * (k - 1)) + (i != 0 ? (N - 1) * (N - 1) : 0) + offset
    end
    offset += 2 * (N - 1) * (N - 1)
    if (jbdy) # On j-normal face
      return (i - 1) + ((N - 1) * (k - 1)) + (j != 0 ? (N - 1) * (N - 1) : 0) + offset
    end
    offset += 2 * (N - 1) * (N - 1)
    # kbdy, On k-normal face
    return (i - 1) + ((N - 1) * (j - 1)) + (k != 0 ? (N - 1) * (N - 1) : 0) + offset
  end

  # nbdy == 0: Body DOF
  offset += 2 * ( (N - 1) * (N - 1) +
                 (N - 1) * (N - 1) +
                 (N - 1) * (N - 1))
  return offset + (i - 1) + (N - 1) * ((j - 1) + (N - 1) * ((k - 1)))
end




let
  N = 5
  nelem = 4
  T = Float64
  dim = 3

  io = open("data$(dim)d.vtu", "w")
  Nq = N+1
  (r, _) = GaussQuadrature.legendre(T, Nq, GaussQuadrature.both)
  s = r
  t = dim == 2 ? [0] : r
  Nqk = length(t)
  Np = Nq * Nq * Nqk


  x = Array{T, 4}(undef, Nq, Nq, Nqk, nelem)
  y = Array{T, 4}(undef, Nq, Nq, Nqk, nelem)
  z = Array{T, 4}(undef, Nq, Nq, Nqk, nelem)

  for e = 1:nelem, k = 1:Nqk, j = 1:Nq, i = 1:Nq
    xoffset = 2e - 1 - nelem
    x[i, j, k, e], y[i, j, k, e], z[i, j, k, e] = r[i]-xoffset, s[j], t[k]
  end
  if dim == 2
    x, y = x + sin.(π * y) / 5, y + exp.(-x.^2)
  else
    x, y, z = x + sin.(π * y) / 5, y + exp.(-hypot.(x, z).^2), z + sin.(π * x) / 5
  end
  d = exp.(sin.(hypot.(x, y, x)))

  write(io,
        """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
          <UnstructuredGrid>
            <Piece NumberOfPoints="$(Np * nelem)" NumberOfCells="$nelem">
              <PointData>
                <DataArray type="Float32" Name="element_number" format="ascii">
        """
       )
  for n = 1:length(d)
    @printf(io, "           %16e\n", div(n-1, Np))
  end
  write(io,
        """
                </DataArray>
                <DataArray type="Float32" Name="data" format="ascii">
        """
       )
  for n = 1:length(d)
    @printf(io, "           %16e\n", d[n])
  end
  write(io,
        """
                </DataArray>
              </PointData>
              <CellData>
              </CellData>
              <Points>
                <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">
        """
       )
  # Fill coordinates
  for n = 1:length(x)
    @printf(io, "           %16e %16e %16e\n", x[n], y[n], z[n])
  end

  write(io,
        """
                </DataArray>
              </Points>
              <Cells>
                <DataArray type="Int64" Name="connectivity" format="ascii">
        """
       )
  # Fill coordinates
  L = LinearIndices((1:Nq, 1:Nq, 1:Nqk))
  for e = 1:nelem
    o = (e-1) * Np - 1
    if dim == 2
      # corners
      @printf(io, "           %d %d %d %d\n",
              o + L[1, 1, 1],
              o + L[end, 1, 1],
              o + L[end, end, 1],
              o + L[1, end, 1])
      # edges
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[n, 1, 1])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[end, n, 1])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[n, end, 1])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[1, n, 1])
      end
      @printf(io, "\n")
      # interior
      for j = 2:Nq-1
        @printf(io, "          ")
        for i = 2:Nq-1
          @printf(io, " %d", o + L[i, j, 1])
        end
        j != Nq-1 && @printf(io, "\n")
      end
      @printf(io, "\n")
    elseif dim == 3
      # corners
      @printf(io, "           %d %d %d %d %d %d %d %d\n",
              o + L[1, 1, 1],
              o + L[end, 1, 1],
              o + L[end, end, 1],
              o + L[1, end, 1],
              o + L[1, 1, end],
              o + L[end, 1, end],
              o + L[end, end, end],
              o + L[1, end, end])
      tmp = 0
      @assert tmp == PointIndexFromIJK(0, 0, 0, N); tmp += 1
      @assert tmp == PointIndexFromIJK(N, 0, 0, N); tmp += 1
      @assert tmp == PointIndexFromIJK(N, N, 0, N); tmp += 1
      @assert tmp == PointIndexFromIJK(0, N, 0, N); tmp += 1
      @assert tmp == PointIndexFromIJK(0, 0, N, N); tmp += 1
      @assert tmp == PointIndexFromIJK(N, 0, N, N); tmp += 1
      @assert tmp == PointIndexFromIJK(N, N, N, N); tmp += 1
      @assert tmp == PointIndexFromIJK(0, N, N, N); tmp += 1
      # edges
      @printf(io, "          ")
      j = 1
      k = 1
      for i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = Nq
      k = 1
      for j = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      j = Nq
      k = 1
      for i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = 1
      k = 1
      for j = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")

      @printf(io, "          ")
      j = 1
      k = Nq
      for i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = Nq
      k = Nq
      for j = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      j = Nq
      k = Nq
      for i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = 1
      k = Nq
      for j = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")

      @printf(io, "          ")
      i = 1
      j = 1
      for k = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = Nq
      j = 1
      for k = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = 1
      j = Nq
      for k = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = Nq
      j = Nq
      for k = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      # faces
      @printf(io, "          ")
      i = 1
      for k = 2:Nq-1, j = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      i = Nq
      for k = 2:Nq-1, j = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      j = 1
      for k = 2:Nq-1, i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      j = Nq
      for k = 2:Nq-1, i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      k = 1
      for j = 2:Nq-1, i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      @printf(io, "          ")
      k = Nq
      for j = 2:Nq-1, i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
      # interior
      @printf(io, "          ")
      for k = 2:Nq-1, j = 2:Nq-1, i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
        @assert tmp == PointIndexFromIJK(i-1, j-1, k-1, N); tmp += 1
      end
      @printf(io, "\n")
    else
      error("invalid dim")
    end
  end
  write(io,
        """
                </DataArray>
                <DataArray type="Int64" Name="offsets" format="ascii">
        """
       )
  for e = 1:nelem
    @printf(io, "           %d\n", Np * e)
  end
  write(io,
        """
                </DataArray>
                <DataArray type="UInt8" Name="types" format="ascii">
        """
       )
  if dim == 2
    for e = 1:nelem
      @printf(io, "           70\n")
    end
  elseif dim == 3
    for e = 1:nelem
      @printf(io, "           72\n")
    end
  else
    error("invalid dim")
  end
  write(io,
        """
                </DataArray>
              </Cells>
            </Piece>
          </UnstructuredGrid>
        </VTKFile>
        """
       )



  close(io)
end
