function ke=kw51_beamMatrix(b,e)
    [ok,ix]=ismember(e(5:7),b.Nodes(:,1)); assert(all(ok));
    p=b.Nodes(ix(1),2:4); q=b.Nodes(ix(2),2:4); ref=b.Nodes(ix(3),2:4);
    L=norm(q-p); ex=(q-p)/L; ez=cross(ex,ref-p); ez=ez/norm(ez); ey=cross(ez,ex);
    rot=[ex;ey;ez]; Q=blkdiag(rot,rot,rot,rot);
    section=b.Sections(b.Sections(:,1)==e(3),2:end);
    material=b.Materials(b.Materials(:,1)==e(4),2:end);
    A=section(1); ky=section(2); kz=section(3); J=section(4); Iy=section(5); Iz=section(6);
    E=material(1); nu=material(2); G=E/(2*(1+nu)); kl=zeros(12);
    kl([1 7],[1 7])=E*A/L*[1 -1;-1 1];
    kl([4 10],[4 10])=G*J/L*[1 -1;-1 1];
    Is=[Iz Iy]; ks=[ky kz]; signs=[1 -1]; indices={[2 6 8 12],[3 5 9 11]};
    for j=1:2
        phi=12*E*Is(j)/(G*A*L^2*ks(j));
        aa=12*E*Is(j)/(L^3*(1+phi)); bb=signs(j)*aa*L/2;
        cc=(4+phi)*E*Is(j)/(L*(1+phi)); dd=(2-phi)*E*Is(j)/(L*(1+phi));
        ii=indices{j}; kl(ii,ii)=[aa bb -aa bb;bb cc -bb dd;-aa -bb aa -bb;bb dd -bb cc];
    end
    ke=Q.'*kl*Q;
end
