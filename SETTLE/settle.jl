export 
    SETTLE,
    apply_constraint!



struct SETTLE{D, B, E}
    dists::D
    is::B
    js::B
    tolerance::E
end

SETTLE(dists, is, js, tolerance=1e-10u"nm") = SETTLE{typeof(dists), typeof(tolerance)}(dists, is, js, tolerance)

function apply_constraints!()
    


end





