using Printf
try
  GaussQuadrature.neither
catch
  include("quadrature.jl")
end

let
  N = 3
  nelem = 4
  T = Float64
  dim = 3

  io = open("proposed_data$(dim)d.vtu", "w")
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
            <Piece NumberOfPoints="$(Nq + Np * nelem)" NumberOfCells="$nelem">
              <PointData>
                <DataArray type="Float32" Name="data" format="ascii">
        """
       )
  for n = 1:Nq
    @printf(io, "           %16e\n", 0)
  end
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
  for n = 1:Nq
    @printf(io, "           %16e %16e %16e\n", r[n], 0, 0)
  end
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
    o = Nq + (e-1) * Np - 1
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

      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[n, 1, end])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[end, n, end])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[n, end, end])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[1, n, end])
      end
      @printf(io, "\n")

      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[1, 1, n])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[end, 1, n])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[1, end, n])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for n = 2:Nq-1
        @printf(io, " %d", o + L[end, end, n])
      end
      @printf(io, "\n")
      # faces
      @printf(io, "          ")
      for m = 2:Nq-1, n = 2:Nq-1
        @printf(io, " %d", o + L[1, n, m])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for m = 2:Nq-1, n = 2:Nq-1
        @printf(io, " %d", o + L[end, n, m])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for m = 2:Nq-1, n = 2:Nq-1
        @printf(io, " %d", o + L[n, 1, m])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for m = 2:Nq-1, n = 2:Nq-1
        @printf(io, " %d", o + L[n, end, m])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for m = 2:Nq-1, n = 2:Nq-1
        @printf(io, " %d", o + L[n, m, 1])
      end
      @printf(io, "\n")
      @printf(io, "          ")
      for m = 2:Nq-1, n = 2:Nq-1
        @printf(io, " %d", o + L[n, m, end])
      end
      @printf(io, "\n")
      # interior
      @printf(io, "          ")
      for k = 2:Nq-1, j = 2:Nq-1, i = 2:Nq-1
        @printf(io, " %d", o + L[i, j, k])
      end
      @printf(io, "\n")
    else
      error("invalid dim")
    end
    # DOF grid
    @printf(io, "          ")
    for n = 0:N
      @printf(io, " %d", n)
    end
    @printf(io, "\n")
  end
  write(io,
        """
                </DataArray>
                <DataArray type="Int64" Name="offsets" format="ascii">
        """
       )
  for e = 1:nelem
    @printf(io, "           %d\n", (Np + Nq) * e)
  end
  write(io,
        """
                </DataArray>
                <DataArray type="UInt8" Name="types" format="ascii">
        """
       )
  if dim == 2
    for e = 1:nelem
      @printf(io, "           75\n")
    end
  elseif dim == 3
    for e = 1:nelem
      @printf(io, "           76\n")
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
