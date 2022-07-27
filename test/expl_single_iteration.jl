#using Molly

#include("../explicit/explicit.jl")

using ConstrainedMolly
using Molly

e = explicit(0.2, 10000, 0.001)

n_atoms = 4

box_size = SVector(2.0, 2.0, 2.0)

coords = place_atoms(n_atoms ÷ 2, box_size, 0.3)
for i in 1:length(coords)
    push!(coords, coords[i] .+ [0.1, 0.0, 0.0])
end


bonds = InteractionList2Atoms(
    collect(1:(n_atoms ÷ 2)),           # First atom indices
    collect((1 + n_atoms ÷ 2):n_atoms), # Second atom indices
    repeat([""], n_atoms ÷ 2),          # Bond types
    [HarmonicBond(b0=0.1u"nm", kb=300_000.0u"kJ * mol^-1 * nm^-2") for i in 1:(n_atoms ÷ 2)],
)

masses = [10 for i in 1:n_atoms]

println("Initial : ")

println(coords)

dt = 0.01

v = [[0.1, 0.1, 0.1] for i in 1:n_atoms]

new_coords = coords + (v .* dt) 

apply_constraint!(bonds, coords, new_coords, masses, dt, e)


