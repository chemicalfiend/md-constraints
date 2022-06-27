using ConstrainedMolly
#using GLMakie


function system_test()
    
    n_atoms = 100
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    
    box_size = SVector(2.0, 2.0, 2.0)u"nm"
    coords = place_atoms(n_atoms, box_size, 0.3u"nm")
    
    temp = 100.0u"K"

    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
    pairwise_inters = (LennardJones(),)
    
    constraint_list = (NoConstraint(),)
    
    sys = ConstrainedSystem(atoms=atoms,
                            pairwise_inters=pairwise_inters,
                            constraints=constraint_list,
                            coords=coords,
                            velocities=velocities,
                            box_size=box_size,
                            loggers=(coords = CoordinateLogger(10),)
                           )
    return sys
end


function ExplicitTest()
    n_atoms = 20
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    
    box_size = SVector(2.0, 2.0, 2.0)u"nm"
    coords = place_atoms(n_atoms ÷ 2, box_size, 0.3u"nm")
    
    for i in 1:length(coords)
        push!(coords, coords[i] .+ [0.1, 0.0, 0.0]u"nm")
    end

    bonds = InteractionList2Atoms(
         collect(1:(n_atoms÷2)),
         collect((1+n_atoms÷2):n_atoms),
         repeat([""], n_atoms÷2),
         [HarmonicBond(b0=0.1u"nm", kb=300_000.0u"kJ * mol^-1 * nm^-2") for i in 1:(n_atoms ÷ 2)],
        )

    specific_inter_lists = (bonds,)

    temp = 100.0u"K"

    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
    pairwise_inters = (LennardJones(),)
    
    exp = explicit(0.01u"nm", bonds, 1e6, 1e-4u"nm")
    constraint_list = (exp,)
    
    sys = ConstrainedSystem(atoms=atoms,
                            pairwise_inters=pairwise_inters,
                            constraints=constraint_list,
                            coords=coords,
                            velocities=velocities,
                            box_size=box_size,
                            loggers=(coords = CoordinateLogger(10),)
                           )
    return sys

end

function ConstrainedVerletTest(sys)
    temp = 100.0u"K"
    box_size = SVector(2.0, 2.0, 2.0)u"nm"
    simulator = ConstrainedVelocityVerlet(dt=0.01u"ps", coupling=AndersenThermostat(temp, 1.0u"ps"),)
    simulate!(sys, simulator, 1_000)
#    visualize(sys.sys.loggers.coords, box_size, "sim_lj.mp4")
end



#system_test()
#ConstrainedVerletTest(system_test())

ConstrainedVerletTest(ExplicitTest())



