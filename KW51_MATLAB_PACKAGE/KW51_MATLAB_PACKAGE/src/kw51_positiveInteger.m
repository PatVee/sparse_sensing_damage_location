function yes=kw51_positiveInteger(x)
    yes=isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x==floor(x);
end
