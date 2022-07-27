using Molly
using ConstrainedMolly
using GLMakie

begin

n_atoms = 20

box_size = SVector(2.0, 2.0, 2.0)

atom_mass = 10.0u"u"

atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]

coords = place_atoms(n_atoms ÷ 2, box_size, 0.3)
for i in 1:length(coords)
    push!(coords, coords[i] .+ [0.1, 0.0, 0.0])
end

bonds = InteractionList2Atoms(
    collect(1:(n_atoms ÷ 2)),
    collect((1 + n_atoms ÷ 2):n_atoms),
    repeat([""], n_atoms ÷ 2),
    [HarmonicBond(b0=0.1u"nm", kb=300_000.0u"kJ * mol^-1 * nm^-2") for i in 1:(n_atoms ÷ 2)],
)

pairwise_inters = (LennardJones(nl_only=true),)

masses = [0.01 for i in 1_atoms]

velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]

e = explicit(0.01, bonds, 1000, 0.001)

constraint_list = (e,)

sys = ConstrainedSystem(atoms=atoms, bonds=bonds, constraints=constraint_list, velocities=velocities, box_size=box_size, coords=coords)

simulator = ConstrainedVelocityVerlet(dt=0.01)

simulate!(sys, simulator, 10)
end
