    nnod=4;
nelem = 1;
ntype=1;
nfunc=1;
type=2;
npe=4;
ndfn=3;
fi=0;
ncon=[1 2 4 3];
xy=[0,0
    5,0
    0,5
    5,5];
x=xy(:,1);
y=xy(:,2);

nebc=9;
iebc=[2,3,4,6,7,8,10,11,12];
vebc=[0,0,0,0,0,0,0,0,0];

E=1*10^7;  NU=0.25; h=0.1;
nnbc=0;

fc=1;
C=E/(1-NU^2)*[1,NU,0;NU,1,0;0,0,(1-NU)/2];   
K=5/6; 
A=C(3,3)*h;
D=(E*h^3)/(12*(1-NU^2)); 




