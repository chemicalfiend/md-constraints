# Bond constraint algorithm implemented with explicit computation of inverse Jacobian at each iteration for finding the Lagrange multipliers of the constrained optimisation.


using LinearAlgebra

export
    explicit,
    apply_constraint!,
    constraint


struct explicit{T, I, E} <: Constraint
    d::T
    maxiter::I
    tol::E
end


#function apply_constraint!(bond_list::InteractionList2Atoms, coords, masses, sim, constraint::explicit, cl::CoordinateLogger)

function apply_constraint!(bond_list::InteractionList2Atoms, coords, masses, dt, constr::explicit)

    # Is coordinate logger needed here?

    σ = []
    m = length(bond_list.is) # Number of constraints = number of bonds
    J = zeros(m, m)  # Initialising the Jacobian

    for i in 1:m
        push!(σ, constraint(coords[bond_list.is[i]], coords[bond_list.js[i]], constr.d) )
    end
    
    for l in 1:constr.maxiter
        coord_i = 0
        coord_j = 0
        for u in 1:m
            for v in 1:m
                J[u, v] = 2 * norm(coords[bond_list.is[u]] - coords[bond_list.js[u]]) * norm( (grad_constraint(coords[bond_list.is[v]], coords[bond_list.js[v]])/masses[v]) - (grad_constraint(coords[bond_list.js[v]], coords[bond_list.is[v]])/masses[v])  )  # TODO : Check the Jacobian !!
            end
        end
        
        println("\n\n\n\n")

        println(J)

        println(σ)

        λ = - J \ σ      
        
        for r in 1:m  # Replace this with enumerate??
            coord_i = coords[bond_list.is[r]]
            coord_j = coords[bond_list.js[r]]
            
            mi = masses[bond_list.is[r]]
            mj = masses[bond_list.js[r]]

            sumi = zeros(3)
            sumj = zeros(3)
            for k in 1:m
                sumi += (dt)^2 * λ[k] * (grad_constraint(coord_i, coord_j)/mi) 
                sumj += (dt)^2 * λ[k] * (grad_constraint(coord_j, coord_i)/mj)
            end
            
            coord_i += sumi
            coord_j += sumj
                        
            for k in 1:m
                push!(σnew, constraint(coord_i, coord_j))
            end
        end


        if (norm(σ - σnew) <= constraint.tol)           
            coords[bond_list.is[r]] = coord_i
            coords[bond_list.js[r]] = coord_j
            return σnew
        end
        
        println("\n\n\n")

        println("Iteration $l")
        println(coords)
        σ = σnew
      
    end

    println("Maximum Iterations Exceeded, σ did not converge")
        
end

function constraint(coord_i, coord_j, d) 
    return (norm(coord_i - coord_j)^2 - (d^2))
end 

function grad_constraint(coord_i, coord_j)
    return 2 * (coord_i - coord_j)
end
