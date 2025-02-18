

"""
    meshopenbook(angle, nbofpages, h)

Create a mesh of 3D open-book polyhedra characterized by the opening angle `angle`, and the number of pages `nbofpages` >= 2. The polyhedron lies between the planes  z = 0 and z = −1.
The definition of this geometry is in: Chandler-Wilde, S. N. and Spence, E. A. Numerische Mathematik , Vol. 150, No. 2, p. 299-371, 2022. 

The target edge size is `h`.
"""
function meshopenbook(angle, nbofpages, h)
    @assert nbofpages ≥ 2
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("openbook")

    # parameters
    θn = angle/(2*nbofpages - 1)
    η = angle^2/(4*π*(2*nbofpages - 1))
    r1 = cos(θn/2 + η)/cos(θn/2)
    r2 = sin(η)/sin(θn/2)

    # top plate (z = 0)
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h, 1)

    for i in 1:nbofpages
        gmsh.model.geo.addPoint(cos((2*i-2)*θn), sin((2*i-2)*θn), 0.0, h, 9*i+1)
        gmsh.model.geo.addLine(1, 9*i+1, 9*i+1)

        gmsh.model.geo.addPoint(cos((2*i-1)*θn), sin((2*i-1)*θn), 0.0, h, 9*i+2)
        gmsh.model.geo.addLine(9*i+1, 9*i+2, 9*i+2)

        gmsh.model.geo.addLine(9*i+2, 1, 9*i+3)
        gmsh.model.geo.addCurveLoop([9*i+1, 9*i+2, 9*i+3], 9*i+1)
        gmsh.model.geo.addPlaneSurface([9*i+1], 9*i+1)
    end

    # bottom plate (z = -1) and sides
    gmsh.model.geo.addPoint(0.0, 0.0, -1.0, h, 2)
    gmsh.model.geo.addLine(2, 1, 9)

    base_curve_loop = Vector{Int}()
    for i in 1:nbofpages
        if i == 1
            gmsh.model.geo.addPoint(r1, 0.0, -1.0, h, 9*i+4)
            gmsh.model.geo.addLine(2, 9*i+4, 9*i+4)
        else
            gmsh.model.geo.addPoint(cos((2*i-2)*θn - η), sin((2*i-2)*θn - η), -1.0, h, 9*i+4)
            gmsh.model.geo.addLine(9*i-3, 9*i+4, 9*i+4)            
        end
        gmsh.model.geo.addLine(9*i+4, 9*i+1, 9*i+7)
        gmsh.model.geo.addCurveLoop([-(9*i), 9*i+4, 9*i+7, -(9*i+1)], 9*i+2)
        gmsh.model.geo.addPlaneSurface([9*i+2], 9*i+2)

        gmsh.model.geo.addPoint(cos((2*i-1)*θn + η), sin((2*i-1)*θn + η), -1.0, h, 9*i+5)
        gmsh.model.geo.addLine(9*i+4, 9*i+5, 9*i+5)
        gmsh.model.geo.addLine(9*i+5, 9*i+2, 9*i+8)
        gmsh.model.geo.addCurveLoop([-(9*i+7), 9*i+5, 9*i+8, -(9*i+2)], 9*i+3)
        gmsh.model.geo.addPlaneSurface([9*i+3], 9*i+3)

        if i == nbofpages
            gmsh.model.geo.addLine(9*i+5, 2, 9*i+6)
            gmsh.model.geo.addCurveLoop([-(9*i+8), 9*i+6, 9, -(9*i+3)], 9*i+4)
        else
            gmsh.model.geo.addPoint(r2*cos((2*i-0.5)*θn), r2*sin((2*i-0.5)*θn), -1.0, h, 9*i+6)
            gmsh.model.geo.addLine(9*i+5, 9*i+6, 9*i+6)
            gmsh.model.geo.addLine(9*i+6, 1, 9*i+9)
            gmsh.model.geo.addCurveLoop([-(9*i+8), 9*i+6, 9*i+9, -(9*i+3)], 9*i+4)
        end
        gmsh.model.geo.addPlaneSurface([9*i+4], 9*i+4)

        append!(base_curve_loop, [-(9*i+4), -(9*i+5), -(9*i+6)])
    end

    gmsh.model.geo.addCurveLoop(base_curve_loop, 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end



"""
    meshuniformtetrahedron(edgelength, h)

Create a mesh of a uniform tetrahedron with edge length `egdelength`. The four vertices are: (1, 1, 1), (-1, -1, 1), (-1, 1, -1) and (1, -1, -1) (rescaled).

The target edge size is `h`.
"""
function meshuniformtetrahedron(edgelength, h)
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("tetrahedron")

    c = edgelength/(2sqrt(2))
    gmsh.model.geo.addPoint(c, c, c, h, 1)
    gmsh.model.geo.addPoint(-c, -c, c, h, 2)
    gmsh.model.geo.addPoint(-c, c, -c, h, 3)
    gmsh.model.geo.addPoint(c, -c, -c, h, 4)

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(1, 3, 2)
    gmsh.model.geo.addLine(1, 4, 3)
    gmsh.model.geo.addLine(2, 3, 4)
    gmsh.model.geo.addLine(2, 4, 5)
    gmsh.model.geo.addLine(3, 4, 6)

    gmsh.model.geo.addCurveLoop([-1, -4, 2], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addCurveLoop([-2, -6, 3], 2)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addCurveLoop([-3, 5, 1], 3)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addCurveLoop([4, -5, 6], 4)
    gmsh.model.geo.addPlaneSurface([4], 4)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end

