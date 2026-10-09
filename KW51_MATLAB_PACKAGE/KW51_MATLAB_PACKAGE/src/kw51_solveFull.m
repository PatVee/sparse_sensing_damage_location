function U=kw51_solveFull(K,scale,factor,F)
    U=scale.*(factor\(scale.*F));
    for j=1:2
        U=U+scale.*(factor\(scale.*(F-K*U)));
    end
end
