using ConstrainedMolly
using GLMakie


function ShakeTest()
    n_atoms = 100
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")

    coords = place_atoms(n_atoms ÷ 2, boundary, 0.3u"nm")

    for i in 1:length(coords)
        push!(coords, coords[i].+ [0.15, 0.0, 0.0]u"nm")
    end

    temp = 100.0u"K"

    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
    
    bonds = InteractionList2Atoms(
        collect(1:(n_atoms÷2)),
        collect((1+n_atoms÷2):n_atoms),
        repeat([""], n_atoms÷2),
        [HarmonicBond(b0=0.1u"nm", kb=300_000.0u"kJ * mol^-1 * nm^-2") for i in 1:(n_atoms ÷ 2)],
       )
    
    specific_inter_lists = (bonds,)

    nb_matrix = trues(n_atoms, n_atoms)

    for i in (1:n_atoms÷2)
        nb_matrix[i, i + (n_atoms ÷ 2)] = false
        nb_matrix[i + (n_atoms ÷ 2), i] = false
    end

    neighbor_finder = DistanceNeighborFinder(nb_matrix=nb_matrix, n_steps=10, dist_cutoff=1.5u"nm")

    sh = SHAKE(0.1u"nm", bonds)
    
    constraint_list = (sh,)

    sys = ConstrainedSystem(atoms=atoms,
                            pairwise_inters=(LennardJones(nl_only=true),),
                            #specific_inter_lists=specific_inter_lists,
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

    return sys
    
end

function SimulateShake(sys)
    n_atoms = 100
    temp = 100.0u"K"
    boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")
    simulator = ConstrainedVelocityVerlet(dt=0.002u"ps", coupling=AndersenThermostat(temp, 1.0u"ps"),)
    
    simulate!(sys, simulator, 1_000)

    visualize(sys.sys.loggers["coords"],
              boundary,
              "sim_diatomic.mp4";
               connections=[(i, i + (n_atoms ÷ 2)) for i in 1:(n_atoms ÷ 2)],
            )

end

SimulateShake(ShakeTest())
