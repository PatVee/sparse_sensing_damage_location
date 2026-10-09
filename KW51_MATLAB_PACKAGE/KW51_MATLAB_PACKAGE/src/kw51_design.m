function D=kw51_design(m,g,t)
    t(~isfinite(t))=m.tref; z=(t(:)-m.tref)/m.tscale; g=g(:);
    D=[ones(size(g)),g,z,g.*z];
end
