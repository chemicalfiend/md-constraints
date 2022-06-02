# Bond constraint algorithm implemented with explicit computation of inverse Jacobian at each iteration for finding the Lagrange multipliers of the constrained optimisation.


export
    explicit,
    apply_constraint!,
    constraint


struct explicit{T, I, E} <: Constraint
    d::T
    maxiter::I
    tol::E
end


function apply_constraint!(sys, sim, constraint::explicit, cl::CoordinateLogger)
    
    bondlist = [] # TODO => bondlist should be set to have a list of all pairs of bonded atoms after the initial unconstrained step has been taken. 

    σ = [constraint(i, j) for (i, j) in bondlist] # TODO => set constraint vector properly
    
    m = length(bondlist) # Number of constraints = number of bonds

    J = Matrix(I, m, m) # Jacobian
    


    # TODO => write the gradi_sigma function to return the gradient for a constraint function
    
    for l in 1:constraint.maxiter
        
        J = [2 * norm(coord_i - coord_j) * norm((gradi_sigma/mi) - (gradj_sigma/mj))] # TODO => Set the Jacobian matrix correctly.

        λ = - J \ σ

        
        
        for (i, j) in bondlist # TODO => Write this properly after bondlist has been defined correctly.
            sumi = zeros(3)
            sumj = zeros(3)

            for k in 1:m
                sumi += (sim.dt)^2 * λ[k] * (gradi_sigma/mi) 
                sumj += (sim.dt)^2 * λ[k] * (gradj_sigma/mj)
            end
            
            coord_i += sumi
            coord_j += sumj
            
            σnew = [constraint(coord_i, coord_j) for (i, j) in bondlist] # TODO => again, set constraints vector correctly.

            if (norm(σ - σnew) <= constraint.tol)
                
                # TODO : Update simulation positions of all the atoms with the extra constraint step.
            end

            σ = σnew

        end
        
    end

    
    println("Maximum Iterations Exceeded, σ did not converge")
        
end


function constraint(coord_i, coord_j, constraint::explicit)
    
    return (norm(coord_i - coord_j)^2 - (constraint.d)^2)

end 
