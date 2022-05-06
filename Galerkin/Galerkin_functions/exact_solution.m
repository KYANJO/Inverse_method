%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.

%---------------------------------------------------------------------%
function qe = exact_solution(coord,nelem,ngl,intma,npoin,icase,time)

%Initialize
qe=zeros(npoin,1);

%Generate Grid Points
for e=1:nelem
    for i=1:ngl
       x=coord(i,e);
       ip=intma(i,e);
       
       qe(ip)=sin(x)*exp(-time);
           
    end
end

      
