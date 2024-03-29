# Gay-Berne Potential for Polyellipsoid interactions

# Implementation as described in 

export GayBerne


struct GayBerne{C, F, E} <: PairwiseInteraction
    cutoff::C
    force_units::F
    energy_units::E
end


@inline @inbounds function force(inter::GayBerne,
                                    dr,
                                    coord_i,
                                    coord_j,
                                    atom_i,
                                    atom_j,
                                    box_size)
    


    


end

@fastmath function force_divr_nocutoff(::GayBerne) 



end


@inline @inbounds function potential_energy()

end


@fastmath function potential()


end



