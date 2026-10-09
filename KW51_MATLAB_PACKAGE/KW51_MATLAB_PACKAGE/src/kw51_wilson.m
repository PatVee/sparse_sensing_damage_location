function [lo,hi]=kw51_wilson(k,n)
    z=1.959963984540054; p=k/n; den=1+z^2/n;
    centre=(p+z^2/(2*n))/den; half=z*sqrt(p.*(1-p)/n+z^2/(4*n^2))/den;
    lo=max(0,centre-half); hi=min(1,centre+half);
end
