% ===================
% 2-D FEM
% ME7620: ADV FINITE ELEMENT ANALYSIS
% Astrit Imeri - T00249444
% % ===================
% ncon - connectivity array
% ndfn - number of unknowns (d.o.f) per node
% nelem - number of elements 
% nnod - number of nodes
% nnbc, nebc: number of known NBC and EBC
% inbc,iebc: location of known NBC and EBC
% vnbc,vebc: values of known NBC and EBC
%ntype - 1: 2nd order ODE, 2: 4th order ODE
%nfunc - 1: linear 2: quadratic
%type - 1 if triangular element, 2 if quadrilateral
%gk - global coefficient matrix, gf- global right hand side vector
clear all;
clc;
% Call the input file


shl1x1l
% shl2x2l
% shl4x4l
% shl2x2q


if type==1 
%     npe=3;
%     if ntype==2
%          ndfn=3;
%     else
%          ndfn=1;sx_c
%     if nfunc==2
%           npe=6;
%     end
%     end
else 
    npe=4;
    if ntype==2
npe=8
    else
        ndfn=3;
    end
end


% End reading data

% Calculate number of total equations to solve

neq = nnod*ndfn;
g_k = zeros(neq,neq);
g_f = zeros(neq,1);
el_x = zeros(npe,1);
el_y = zeros(npe,1);
sx_c=zeros(nelem,1);
sy_c=zeros(nelem,1);
sxy_c=zeros(nelem,1);
xc_c=zeros(nelem,1);
yc_c=zeros(nelem,1);


% Start loop over number of elements to calculate the element matrices
% and assemble them into the global matrices
%

for n=1:nelem
    
    for i=1:npe
        el_x(i) = x(ncon(n,i));
        el_y(i) = y(ncon(n,i));   
    end
 
    el_x;
    el_y;
    if type==1
    [el_k,el_f] = elkk2d(npe,fc,el_x,el_y,ey,ex);

    else
   [el_k,el_f] = el_kk_quad_shear_lock(ntype,npe,C,fi,fc,D,NU,K,A,el_x,el_y,nfunc,ndfn);
    end
    % Assembly
    for i = 1:npe
        nr = (ncon(n,i)-1)*ndfn;
        for ii = 1:ndfn
  nr = nr + 1;
            l = (i-1)*ndfn+ii;
            g_f(nr) = g_f(nr)+el_f(l);
            for j = 1:npe
                nc = (ncon(n,j)-1)*ndfn;
                for jj = 1:ndfn
                    m = (j-1)*ndfn+jj;
                    nc = nc+1;
                    g_k(nr,nc)=g_k(nr,nc)+el_k(l,m);
                end
            end
        end
    end
end

% Apply known forces

for n = 1:nnbc
     nb = inbc(n);
     g_f(nb) = g_f(nb)+vnbc(n);
 end

% Apply known displacements (1-0) method

for nj = 1:nebc
    j = iebc(nj);
    for k=1:neq
        if k~=j
            g_f(k) = g_f(k)-g_k(k,j)*vebc(nj);
            g_k(k,j) = 0;
            g_k(j,k) = 0;
        else
            g_k(j,j) = 1;
            g_f(j) = vebc(nj);
        end
    end
end
sol = g_k\g_f
%  scatter(x,y)
%     hold on
xx=zeros(nnod,1);yy=xx;
for gnode=1:nnod
    xx(gnode)=x(gnode)+10^2*sol((gnode-1)*ndfn+1);
    yy(gnode)=y(gnode)+10^2*sol((gnode-1)*ndfn+2);
end
% scatter(xx,yy)


sigmax=zeros(nnod,1);
sigmay=zeros(nnod,1);
sigmaxy=zeros(nnod,1);
sw=zeros(nnod,1);
for n=1:nelem
    eta=0.0;
    xi=0.0;
    
    
    for i=1:npe
        nn=ncon(n,i);
        
        el_x(i)=x(nn);
        el_y(i)=y(nn);
        
        ksi_x(i)=sol(ndfn*nn-1);
        ksi_y(i)=sol(ndfn*nn);
        
    end
    
    [sf,dsfx,dsfy,J] = shape_2d_shear_lock(nfunc,npe,xi,eta,el_x,el_y);
    
    sx=0.0;
    sy=0.0;
    sxy=0.0;
    
    xc=0;
    yc=0;
    for i=1:npe
        sx=sx+(C(1,1)*ksi_x(i)*dsfx(i)+C(1,2)*ksi_y(i)*dsfy(i))*h/2;
        sy=sy+(C(1,2)*ksi_x(i)*dsfx(i)+C(2,2)*ksi_y(i)*dsfy(i))*h/2;
        sxy=sxy+(C(3,3)*ksi_x(i)*dsfx(i)+C(3,3)*ksi_y(i)*dsfy(i))*h/2;
        
         
        xc=xc+el_x(i)*sf(i);
        yc=yc+el_y(i)*sf(i);
         
    end
    xc_c(n)=xc;
    yc_c(n)=yc;
    
    sx_c(n)=sx;
    sy_c(n)=sy;
    sxy_c(n)=sxy;
    for k=1:npe
        npk=ncon(n,k);
        xk = el_x(k);
        yk = el_y(k); 
        
        weight=1/(sqrt((xc-xk)^2+(yc-yk)^2));
        sigmax(npk)=sigmax(npk)+weight*sx;
        sigmay(npk)=sigmay(npk)+weight*sy;
        sigmaxy(npk)=sigmaxy(npk)+weight*sxy;
        
        sw(npk)=sw(npk)+weight;
    end
end

for i=1:nnod
    sigmax(i)=sigmax(i)/sw(i);
    sigmay(i)=sigmay(i)/sw(i);
    sigmaxy(i)=sigmaxy(i)/sw(i);
    
    
end
svm=zeros(nnod,1);
for i=1:nnod
    sigma=[sigmax(i), sigmaxy(i);
        sigmaxy(i),sigmay(i)];
    sigmap=eig(sigma);
    if sigmap(1) > sigmap(2)
        svm(i)=sigmap(1);
    else
        svm(i)=sigmap(2);
    end
end
sol1 = sol(1:2:end);

r=sol(1:2:end);
c=sol(2:2:end);

% figure('Renderer', 'painters', 'Position', [10 10 1000 500])
% F=scatteredInterpolant(x,y,r);
% [xgrid,ygrid] = meshgrid(linspace(min(x),max(x)),linspace(min(y),max(y)));
% solgrid=F(xgrid,ygrid);
% contourf(xgrid,ygrid,solgrid)
% colorbar 