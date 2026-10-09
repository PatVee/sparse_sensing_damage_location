function rows=kw51_localRows(e,dof)
    keys=reshape((e(5:6).'*100+(1:6)).',[],1);
    [ok,rows]=ismember(round(keys),dof);
    assert(all(ok),'Missing beam degrees of freedom.');
end
