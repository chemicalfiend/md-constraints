struct SHAKE{D, B, I, E} <: Constraint
    distance::D
    bond_list::B
    maxiter::I
    tol::E
end


function apply_constraint!(sys, coords, new_coords, dt, constraint::SHAKE)
    


end


