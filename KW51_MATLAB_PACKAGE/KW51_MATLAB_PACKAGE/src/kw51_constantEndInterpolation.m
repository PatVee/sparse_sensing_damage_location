function x=kw51_constantEndInterpolation(x,finite)
    ix=find(finite); original=x(ix);
    x=interp1(ix,original,(1:numel(x)).','linear');
    x(1:ix(1)-1)=original(1); x(ix(end)+1:end)=original(end);
end
