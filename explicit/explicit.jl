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


function apply_constraint!(bond_list::InteractionList2Atoms, coords, masses, sim, constraint::explicit, cl::CoordinateLogger)

    # Is coordinate logger needed here?

    σ = []
    m = length(bond_list.is) # Number of constraints = number of bonds
    J = Matrix(I, m, m) # Initialising the Jacobian

    for i in 1:m
        push!(σ, constraint(coords[bond_list.is[i]], coords[bond_list.js[i]] ) )
    end
    
    for l in 1:constraint.maxiter
        coord_i = 0
        coord_j = 0
        for u in 1:m
            for v in 1:m
                J[u][v] = 2 * norm(coords[bond_list.is[u]] - coords[bond_list.js[u]]) * norm( (grad_constraint(coords[bond_list.is[v], coods[bond_list.js[v])/masses[i]) - (grad_constraint(coords[bond_list.js[v], coods[bond_list.is[v])/masses[j])  )  # TODO : Check the Jacobian !!
            end
        end

        λ = - J \ σ      
        
        for i in 1:m   # TODO : Replace 'i' here with something else and change accordingly
            coord_i = coords[bond_list.is[i]]
            coord_j = coords[bond_list.js[i]]
            
            mi = masses[bond_list.is[i]]
            mj = masses[bond_list.js[i]]

            sumi = zeros(3)
            sumj = zeros(3)
            for k in 1:m
                sumi += (sim.dt)^2 * λ[k] * (grad_constraint(coord_i, coord_j)/mi) 
                sumj += (sim.dt)^2 * λ[k] * (grad_constraint(coord_j, coord_i)/mj)
            end
            
            coord_i += sumi
            coord_j += sumj
                        
            for k in 1:m
                push!(σnew, constraint(coord_i, coord_j))
            end
        end


        if (norm(σ - σnew) <= constraint.tol)           
            coords[bond_list.is[i]] = coord_i
            coords[bond_list.js[i]] = coord_j
            return σnew
        end

        σ = σnew
      
    end

    println("Maximum Iterations Exceeded, σ did not converge")
        
end

function constraint(coord_i, coord_j, constraint::explicit) 
    return (norm(coord_i - coord_j)^2 - (constraint.d)^2)
end 

function grad_constraint(coord_i, coord_j)
    return 2 * (coord_i - coord_j)
end
