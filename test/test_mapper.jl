using CompScienceMeshes
using Test

#include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))

sphere = readmesh(joinpath(dirname(@__FILE__),"assets","sphere2.in"))
index_in_sphere = mapper(sphere)
for (i,_c) in enumerate(cells(sphere))
    @eval @test index_in_sphere[$_c] == $i
end

isdownbelow(c) = center(chart(sphere,c))[3] < 0
southern_hemisphere = submesh(isdownbelow, sphere)

# build the restriction operator
R = spzeros(Int, numcells(southern_hemisphere), numcells(sphere))
for (i,_c) in enumerate(cells(southern_hemisphere))
    j = index_in_sphere[_c]
    R[i,j] = 1
end

@test nnz(R) == numcells(southern_hemisphere)

rect = meshrectangle(1.0, 1.0, 1.0, 3)
@test numcells(rect) == 2
tri = Mesh(rect.vertices, rect.faces[1:1])
@test numcells(tri) == 1

R = restriction(tri, rect)
@test R == [1 0]
