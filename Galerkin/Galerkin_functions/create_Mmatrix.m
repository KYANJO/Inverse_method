%---------------------------------------------------------------------%
%This function computes the global mass matrix for either CG or DG.
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi)

%Initialize
Mmatrix=zeros(npoin,npoin);

for e=1:nelem

    %Store Coordinates
    x = coord(:,e);

    dx=x(ngl)-x(1);
    jac=dx/2;

    %Do LGL Integration
    for l=1:nq
        wq=wnq(l)*jac;
        for i=1:ngl
            ip=intma(i,e);
            h_i=psi(i,l);
            for j=1:ngl
                jp=intma(j,e);
                h_j=psi(j,l);
                Mmatrix(ip,jp)=Mmatrix(ip,jp) + wq*h_i*h_j;
            end %j
        end %i
    end %l
end %e



      