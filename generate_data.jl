using Printf
try
  GaussQuadrature.neither
catch
  include("quadrature.jl")
end

let
  N = 3
  nelem = 1
  dim = 2
  T = Float64

  io = open("data2d.vtu", "w")
  Nq = N+1
  Nqk = dim == 2 ? 1 : Nq
  Np = Nq * Nq * Nqk
  (r, _) = GaussQuadrature.legendre(T, Nq, GaussQuadrature.both)
  s = copy(r)
  t = [0]


  x = Array{T, 4}(undef, Nq, Nq, Nqk, nelem)
  y = Array{T, 4}(undef, Nq, Nq, Nqk, nelem)
  z = Array{T, 4}(undef, Nq, Nq, Nqk, nelem)

  for e = 1:nelem, k = 1:Nqk, j = 1:Nq, i = 1:Nq
    x[i, j, k, e], y[i, j, k, e], z[i, j, k, e] = r[i], s[j], t[k]
  end
  # z = sin.(π * x) .* cos.(π * y) .* cos.(π * z)
  d = exp.(-hypot.(x, y, z)*10)


  write(io,
        """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
          <UnstructuredGrid>
            <Piece NumberOfPoints="$(Np * nelem)" NumberOfCells="$nelem">
              <PointData>
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
    # elseif dim == 3
      # corners
      # edges
      # faces
      # interior
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
