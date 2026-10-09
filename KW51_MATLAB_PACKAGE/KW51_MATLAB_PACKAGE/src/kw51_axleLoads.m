function F=kw51_axleLoads(b,dof,x,track)
    allGroups=sort([2820:10:2950,2835:10:2945]);
    groups=[2830:10:2940,2835:10:2945];
    cells=zeros(numel(groups)*75,5); rr=0;
    for g=groups
        gi=find(allGroups==g);
        for i=2:76
            node=g*100+i;
            ids=[node,allGroups(gi+1)*100+i,allGroups(gi-1)*100+i,node+1,node-1];
            [ok,ix]=ismember(ids,b.Nodes(:,1)); assert(all(ok));
            p=b.Nodes(ix,2:4); rr=rr+1;
            cells(rr,:)=[node,(p(1,1)+p(5,1))/2,(p(1,1)+p(4,1))/2, ...
                (p(1,2)+p(3,2))/2,(p(1,2)+p(2,2))/2];
        end
    end
    if strcmp(track,'B'), yp=[0 1.56 0]+3.24; else, yp=[0 1.57 0]+7.24; end
    yy=polyval(polyfit([0 58.158 115.658],yp,2),x);
    dx=max(0,min(cells(:,3),x+2.4)-max(cells(:,2),x-2.4));
    dy=max(0,min(cells(:,5),yy+1.3175)-max(cells(:,4),yy-1.3175));
    a=dx.*dy;
    extent=max(0,min(min((x+2.4)/4.8,(116.316-(x-2.4))/4.8),1));
    a=a./max(sum(a,1),1e-30).*extent;
    [ok,ixz]=ismember(round(cells(:,1)*100+3),dof); assert(all(ok));
    F=-b.T(ixz,:).'*sparse(a);
end
