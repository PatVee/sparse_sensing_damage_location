function row=kw51_zoneRow(name,family,q,index)
    w=max(q.elementLength_m,1e-9);
    xyz=sum(q{:,{'x_m','y_m','z_m'}}.*w,1)/sum(w);
    row={index,string(name),string(family),xyz(1),xyz(2),xyz(3), ...
        sum(q.elementLength_m),height(q),strjoin(string(q.elementID),';')};
end
