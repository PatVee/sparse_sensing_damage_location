function [zones, passage, audit] = kw51_build_fe_library(feFile, candidateFile, verbose)
    regions = readtable(candidateFile, 'TextType', 'string');
    zones = kw51_defineZones(regions);
    loaded = load(feFile, 'Data');
    b = loaded.Data.Base;
    K = b.Kr; T = b.T;
    dof = round(double(b.DOF(:))*100);
    nd = size(K,1);
    scale = 1./sqrt(diag(K));
    D = spdiags(scale,0,nd,nd);
    factor = decomposition(D*K*D,'lu');
    S = sparse(12,numel(dof));
    auditNames = strings(0,1); auditValues = zeros(0,1);
    for j=1:12
        e=b.Elements(b.Elements(:,1)==b.gaugeHostElementIDs(j),:);
        rows=kw51_localRows(e,dof);
        A=b.Phi(rows,:).'; target=b.strainModeShapes(j,:).';
        coef=pinv(A,1e-13*norm(A,2))*target;
        S(j,rows)=coef.';
        auditNames(end+1,1)="sensor_modal_image_relative_error_"+j;
        auditValues(end+1,1)=norm(A*coef-target)/norm(target);
    end
    Sr=S*T; active=[1 2 3 5 6 7 8 9 10 11 12];
    x=-3:0.5:119.5;
    axles=[0 2.3 18.6 20.9 26.1 28.4 44.7 47 52.2 54.5 70.8 73.1];
    loads=[132 131.4 116.8 116.4 94.7 94.6 94.6 94.7 116.4 116.8 131.4 132]*1000;
    front=linspace(-3,116.316+axles(end)+3,128);
    L=kw51_loadTransform(x,front,axles,loads);
    U=cell(2,1); H0=zeros(2,11,128);
    for g=1:2
        tracks={'A','B'};
        F=kw51_axleLoads(b,dof,x,tracks{g});
        U{g}=kw51_solveFull(K,scale,factor,full(F))*L;
        H0(g,:,:)=Sr(active,:)*U{g};
        residual=norm(K*U{g}-F*L,'fro')/norm(F*L,'fro');
        auditNames(end+1,1)="equilibrium_relative_error_"+tracks{g};
        auditValues(end+1,1)=residual;
    end
    losses=[.1 .2 .3 .4 .5];
    response=zeros(2,height(zones),5,11,128,'single');
    hosts=double(b.gaugeHostElementIDs(active));
    for z=1:height(zones)
        ids=str2double(split(zones.element_ids(z),';'));
        B=[];
        for eID=ids(:).'
            e=b.Elements(b.Elements(:,1)==eID,:);
            ke=kw51_beamMatrix(b,e);
            [V,E]=eig((ke+ke.')/2,'vector');
            ok=E>max(E)*1e-10;
            assert(sum(ok)==6,'Beam update must have six deformational modes.');
            Be=b.T(kw51_localRows(e,dof),:).'*(V(:,ok).*sqrt(E(ok)).');
            B=[B,full(Be)]; %#ok<AGROW>
        end
        Z=kw51_solveFull(K,scale,factor,B);
        C=(B.'*Z+Z.'*B)/2;
        W=Sr(active,:)*Z;
        for g=1:2
            R=B.'*U{g};
            for si=1:5
                loss=losses(si);
                h=squeeze(H0(g,:,:))+W*((eye(size(C))-loss*C)\(loss*R));
                h(ismember(hosts,ids),:)=h(ismember(hosts,ids),:)*(1-loss);
                response(g,z,si,:,:)=single(h);
            end
        end
        if verbose && mod(z,10)==0, fprintf('  FE zones %d/%d\n',z,height(zones)); end
    end
    passage=struct('response',response,'base',H0,'losses',losses,'front',front);
    audit=table(auditNames,auditValues,'VariableNames',{'test','value'});
end
