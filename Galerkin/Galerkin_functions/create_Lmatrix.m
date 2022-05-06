%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.

%---------------------------------------------------------------------%
function L = create_Lmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,dpsi)

%Initialize

L = zeros(npoin,npoin);

for e=1:nelem
   
   %Store Coordinates
   
   x = coord(:,e);
   
   dx = x(ngl)-x(1);
   
   for l = 1:nq
       
      wq = wnq(l);
            
      %Form Derivative

      for j = 1:ngl
          jp = intma(j,e);
          dhdx_j = dpsi(j,l);
      
          %Form RHS
          for i = 1:ngl
             dhdx_i = dpsi(i,l);
             ip = intma(i,e);
             
             L(ip,jp) = L(ip,jp) - 2*wq*dhdx_i*dhdx_j/dx;
             
          end
      end %i
   end %l
   
end %e

