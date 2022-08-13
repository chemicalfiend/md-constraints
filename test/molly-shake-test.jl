using Molly
#using GLMakie

using LinearAlgebra

function ShakeTest()
    n_atoms = 100
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")

    coords = place_atoms(n_atoms ÷ 2, boundary, min_dist=0.3u"nm")

    for i in 1:length(coords)
        push!(coords, coords[i].+ [0.15, 0.0, 0.0]u"nm")
    end

    temp = 100.0u"K"

    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
    #=
    bonds = InteractionList2Atoms(
        collect(1:(n_atoms÷2)),
        collect((1+n_atoms÷2):n_atoms),
        repeat([""], n_atoms÷2),
        [],
       )
    =#

    nb_matrix = trues(n_atoms, n_atoms)

    for i in (1:n_atoms÷2)
        nb_matrix[i, i + (n_atoms ÷ 2)] = false
        nb_matrix[i + (n_atoms ÷ 2), i] = false
    end

    neighbor_finder = DistanceNeighborFinder(nb_matrix=nb_matrix, n_steps=10, dist_cutoff=1.5u"nm")
    
    bond_lengths = [0.1u"nm" for i in 1:(n_atoms÷2)]

    sh = SHAKE(bond_lengths, collect(1:(n_atoms÷2)), collect((1+n_atoms÷2):n_atoms))
    
    constraint_list = (sh,)

    sys = System(atoms=atoms,
                 pairwise_inters=(LennardJones(nl_only=true),),

                            constraints=constraint_list,
                            coords=coords,
                            velocities=velocities,
                            boundary=boundary,
                            neighbor_finder=neighbor_finder,
                            loggers=Dict(
                                         "temp" => TemperatureLogger(10),
                                         "coords" => CoordinateLogger(10),
                                        )
                           )
    #n_atoms = 100
    #temp = 100.0u"K"
    #boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")
    
    for i in 1:length(sys.coords)
        sys.coords[i] += [rand()*0.01, rand()*0.01, rand()*0.01]u"nm"
    end
    
    old_coords = sys.coords

    apply_constraint!(sys, sh, old_coords, 0.002u"ps")

    lengths = []
    for r in 1:length(sh.is)
        push!(lengths, norm(vector(sys.coords[sh.is[r]], sys.coords[sh.js[r]], sys.boundary)))
    end
    
    println("TEST 1 : \n")

    print(lengths.-0.1u"nm")
    
    if abs(maximum(lengths.-0.1u"nm")) < 1e-7u"nm"
        println("\n SUCCESS!")
    else
        println("\n FAILED :(")
    end


    simulator = VelocityVerlet(dt=0.002u"ps", coupling=AndersenThermostat(temp, 1.0u"ps"),)

    simulate!(sys, simulator, 1_000)
    
    lengths = []
    for r in 1:length(sh.is)
        push!(lengths, norm(vector(sys.coords[sh.is[r]], sys.coords[sh.js[r]], sys.boundary)))
    end
    
    println("TEST 2 \n")

    print(lengths.-0.1u"nm")
    
    if abs(maximum(lengths.-0.1u"nm")) < 1e-7u"nm"
        println("\n SUCCESS!")
    else
        println("\n FAILED :(")
    end


    #=
    visualize(sys.loggers["coords"],
              boundary,
              "sim_diatomic.mp4";
               connections=[(i, i + (n_atoms ÷ 2)) for i in 1:(n_atoms ÷ 2)],
               connection_frames = values(sys.loggers["bonds"])
            )
    =#
end

ShakeTest()
