# Explicit calculation of matrix inverse for constrained motion of atoms. Refer Ryckaert, Cicottti, Berendsen (1977) Section 2 and 3.


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
#        for u in 1:m
 #           for v in 1:m
  #              J[u, v] = 2 * norm(coords[bond_list.is[v]] - coords[bond_list.js[v]]) * norm( (grad_constraint(coords[bond_list.is[u]], coords[bond_list.js[u]])/masses[u]) - (grad_constraint(coords[bond_list.js[u]], coords[bond_list.is[u]])/masses[u])  )  # TODO : Check the Jacobian !!

#            end
 #       end

        J[1, 1] = 2 * norm(coords[bond_list.is[1]] - coords[bond_list.js[1]]) * norm((grad_constraint(coords[bond_list.is[1]], coords[bond_list.js[1]])/masses[1]) - (grad_constraint(coords[bond_list.js[1]], coords[bond_list.is[1]])/masses[1]))
        

        J[2, 2] = 2 * norm(coords[bond_list.is[2]] - coords[bond_list.js[2]]) * norm((grad_constraint(coords[bond_list.is[2]], coords[bond_list.js[1]])/masses[2]) - (grad_constraint(coords[bond_list.js[2]], coords[bond_list.is[2]])/masses[2])  )


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
            
            σnew = []

            for k in 1:m
                push!(σnew, constraint(coord_i, coord_j, constr.d))
            end

            if (norm(σ - σnew) <= constr.tol)
                coords[bond_list.is[r]] = coord_i
                coords[bond_list.js[r]] = coord_j
                return σnew
            end
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
