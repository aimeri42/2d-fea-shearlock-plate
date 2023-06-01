function [el_k,el_f] = el_kk_quad_shear_lock(ntype,npe,C,fi,fc,D,NU,K,A,el_x,el_y,nfunc,ndfn)

% Gauss points array
gauss = [0.0,-1/sqrt(3),-sqrt(3/5),-0.86113631;
    0.0,1/sqrt(3),0.0,-0.33998104;
    0.0,0.0,sqrt(3/5),0.33998104;
    0.0,0.0,0.0,0.86113631];
% Gauss weights array

weight = [2.0,1.0,5/9,0.34785485;
    0.0,1.0,8/9,0.65214515;...
    0.0,0.0,5/9,0.65214515;
    0.0,0.0,0.0,0.34785485];
% Gauss points array for quad element 4 nodes

if nfunc==2
    ngp=3;
    npe=8;
     npea=npe*ndfn;
else
    ngp=2;
    
    npea=npe*ndfn;
end


el_f = zeros(npea,1);
el_k = zeros(npea,npea);
el_k_11=zeros(npe,npe);
el_k_12=zeros(npe,npe);
el_k_13=zeros(npe,npe);
el_k_21=zeros(npe,npe);
el_k_22=zeros(npe,npe);
el_k_23=zeros(npe,npe);
el_k_31=zeros(npe,npe);
el_k_32=zeros(npe,npe);
el_k_33=zeros(npe,npe);
for ngi =1:ngp
    for ngj =1:ngp
        xi = gauss(ngi,ngp);
        eta = gauss(ngj,ngp);
        
        [sf,dsfx,dsfy,J] = shape_2d_shear_lock(nfunc,npe,xi,eta,el_x,el_y);
        %%
        for i=1:npe
            el_f((i-1)*ndfn+1) = el_f((i-1)*ndfn+1)+fc*sf(i)*weight(ngi,ngp)*J*weight(ngj,ngp);
        end
        
        if ntype==1
            for i=1:npe
                for j=1:npe
                    el_k_22(i,j) = el_k_22(i,j)+(D*dsfx(i)*dsfx(j)...
                        +0.5*(1-NU)*D*dsfy(i)*dsfy(j))*J*weight(ngi,ngp)*weight(ngj,ngp);
                    el_k_23(i,j)=el_k_23(i,j)+(NU*D*dsfx(i)*dsfy(j)+...
                        0.5*(1-NU)*D*dsfy(i)*dsfx(j))*J*weight(ngi,ngp)*weight(ngj,ngp);
                    el_k_32(i,j)=el_k_32(i,j)+(NU*D*dsfy(i)*dsfx(j)...
                        +0.5*(1-NU)*D*dsfx(i)*dsfy(j))*J*weight(ngi,ngp)*weight(ngj,ngp);
                    el_k_33(i,j)=el_k_33(i,j)+(0.5*(1-NU)*D*dsfx(i)*dsfx(j)...
                        +D*dsfy(i)*dsfy(j))*J*weight(ngi,ngp)*weight(ngj,ngp);
                    
                    
                end
            end
        end
    end
end

if fi==1
ngpr=ngp;
else
    ngpr=ngp-1;
end

for ngi=1:ngpr
    for ngj=1:ngpr
        
        xi = gauss(ngi,ngpr);
        eta = gauss(ngj,ngpr);
        
        [sf,dsfx,dsfy,J] = shape_2d_shear_lock(nfunc,npe,xi,eta,el_x,el_y);
        if ntype==1
            for i=1:npe
                for j=1:npe
                    
                    el_k_11(i,j)=el_k_11(i,j)+(K*A*dsfx(i)*dsfx(j)+...
                        K*A*dsfy(i)*dsfy(j))*J*weight(ngi,ngpr)*weight(ngj,ngpr);
                    el_k_12(i,j)=el_k_12(i,j)+(K*A*dsfx(i)*sf(j))...
                        *J*weight(ngi,ngpr)*weight(ngj,ngpr);
                    el_k_13(i,j)=el_k_13(i,j)+(K*A*dsfy(i)*sf(j))...
                        *J*weight(ngi,ngpr)*weight(ngj,ngpr);
                    el_k_21(i,j)=el_k_21(i,j)+K*A*sf(i)*dsfx(j)...
                        *J*weight(ngi,ngpr)*weight(ngj,ngpr);
                    el_k_22(i,j)=el_k_22(i,j)+K*A*sf(i)*sf(j)...
                        *J*weight(ngi,ngpr)*weight(ngj,ngpr);
                    el_k_31(i,j)=el_k_31(i,j)+K*A*sf(i)*dsfy(j)...
                        *J*weight(ngi,ngpr)*weight(ngj,ngpr);
                    el_k_33(i,j)=el_k_33(i,j)+K*A*sf(i)*sf(j)...
                        *J*weight(ngi,ngpr)*weight(ngj,ngpr);
                end
            end
        end
        
    end
end

% el_k=[el_k_11,el_k_12,el_k_13;
%     el_k_21,el_k_22,el_k_23;
%     el_k_31,el_k_32,el_k_33];
ii=1;
for i=1:npe
    jj=1;
    for j=1:npe
        el_k(ii,jj)=el_k_11(i,j);
        el_k(ii,jj+1)=el_k_12(i,j);
        el_k(ii,jj+2)=el_k_13(i,j);
        el_k(ii+1,jj)=el_k_21(i,j);
        el_k(ii+1,jj+1)=el_k_22(i,j);
        el_k(ii+1,jj+2)=el_k_23(i,j);
        el_k(ii+2,jj)=el_k_31(i,j);
        el_k(ii+2,jj+1)=el_k_32(i,j);
        el_k(ii+2,jj+2)=el_k_33(i,j);
        
        jj=ndfn*j+1;
    end
    ii=ndfn*i+1;
end

end