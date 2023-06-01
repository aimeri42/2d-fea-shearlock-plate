function [sf,dsfx,dsfy,J]=shape_2d_shear_lock(nfunc,npe,xi,eta,el_x,el_y)
%
sf = zeros(npe,1);
dsf = zeros(npe,1);
ddsf = zeros(npe,1);
gdsf = zeros(npe,1);
gddsf = zeros(npe,1);
dsfxi = zeros(npe,1);
dsfeta = zeros(npe,1);

   
if nfunc==1

    sf(1) = (1/4)*(1-xi)*(1-eta);
    sf(2) = (1/4)*(1+xi)*(1-eta);
    sf(3) = (1/4)*(1+xi)*(1+eta);
    sf(4) = (1/4)*(1-xi)*(1+eta);

    dsfxi(1)=(1/4)*(eta-1);
    dsfxi(2)=(1/4)*(1-eta);
    dsfxi(3)=(1/4)*(1+eta);
    dsfxi(4)=(1/4)*(-1-eta);

    dsfeta(1)=(1/4)*(xi-1);
    dsfeta(2)=(1/4)*(-1-xi);
    dsfeta(3)=(1/4)*(1+xi);
    dsfeta(4)=(1/4)*(1-xi);

else
 %% Shape functions of 8 node rectangular element
 
    sf(1)=(1/4)*(1-xi)*(1-eta)*(-1-xi-eta);
    sf(2)=(1/4)*(1+xi)*(1-eta)*(-1+xi-eta);
    sf(3)=(1/4)*(1+xi)*(1+eta)*(-1+xi+eta);
    sf(4)=(1/4)*(1-xi)*(1+eta)*(-1-xi+eta);
    sf(5)=(1/2)*(1-xi^2)*(1-eta);
    sf(6)=(1/2)*(1+xi)*(1-eta^2);
    sf(7)=(1/2)*(1-xi^2)*(1+eta);
    sf(8)=(1/2)*(1-xi)*(1-eta^2);
    
    
    dsfxi(1)=-(1/4)*(1-eta)*(-1-xi-eta)-(1/4)*(1-xi)*(1-eta);
    dsfxi(2)=(1/4)*(1-eta)*(-1+xi-eta)+(1/4)*(1+xi)*(1-eta);
    dsfxi(3)=(1/4)*(1+eta)*(-1+xi+eta)+(1/4)*(1+xi)*(1+eta);
    dsfxi(4)=-(1/4)*(1+eta)*(-1-xi+eta)-(1/4)*(1-xi)*(1+eta);
    dsfxi(5)=-xi*(1-eta);
    dsfxi(6)= -(1/2)*eta^2+1/2;
    dsfxi(7)=-xi*(1+eta);
    dsfxi(8)=(1/2)*eta^2-1/2;
    
    dsfeta(1)=-(1/4)*(1-xi)*(-1-xi-eta)-(1/4)*(1-xi)*(1-eta);
    dsfeta(2)=-(1/4)*(1+xi)*(-1+xi-eta)-(1/4)*(1+xi)*(1-eta);
    dsfeta(3)=(1/4)*(1+xi)*(-1+xi+eta)+(1/4)*(1+xi)*(1+eta);
    dsfeta(4)=(1/4)*(-1-xi+eta)*(1-xi)+(1/4)*(1-xi)*(1+eta);
    dsfeta(5)=(1/2)*xi^2-1/2;
    dsfeta(6)=-(1+xi)*eta;
    dsfeta(7)=-(1/2)*xi^2+1/2;
    dsfeta(8)= -(1-xi)*eta;
    
end

% Calculate Jacobian matrix and its determinant

rj=zeros(2,2);

for i=1:npe
    
    rj(1,1)= rj(1,1) + dsfxi(i)*el_x(i);
    rj(1,2)= rj(1,2) + dsfxi(i)*el_y(i);
    rj(2,1)= rj(2,1) + dsfeta(i)*el_x(i);
    rj(2,2)= rj(2,2) + dsfeta(i)*el_y(i);
end   
J=det(rj);
rjinv = inv(rj);

for i=1:npe
    dsfx(i) = rjinv(1,1)*dsfxi(i)+rjinv(1,2)*dsfeta(i);
    dsfy(i) = rjinv(2,1)*dsfxi(i)+rjinv(2,2)*dsfeta(i);
end

end

