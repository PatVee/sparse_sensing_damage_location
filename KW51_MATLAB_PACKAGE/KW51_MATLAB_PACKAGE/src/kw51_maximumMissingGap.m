function n=kw51_maximumMissingGap(finite)
    edges=diff([true;finite(:);true]); starts=find(edges==-1); ends=find(edges==1);
    if isempty(starts), n=0; else, n=max(ends-starts); end
end
