nnod=25;
nelem = 16;
ntype=1;
nfunc=1;
type=2;
npe=4;
ndfn=3;
fi=0;
ncon = [1 2  7 6
    2 3 8 7
    3 4 9 8
    4 5 10 9
    6 7 12 11
    7 8 13 12
    8 9 14 13
    9 10 15 14
    11 12 17 16
    12 13 18 17
    13 14 19 18
    14 15 20 19
    16 17 22 21
    17 18 23 22
    18 19 24 23
    19 20 25 24];

xy=[0,0
1.25,0
2.5,0
3.75,0
5,0
0,1.25
1.25,1.25
2.5,1.25
3.75,1.25
5,1.25
0,2.5
1.25,2.5
2.5,2.5
3.75,2.5
5,2.5
0,3.75
1.25,3.75
2.5,3.75
3.75,3.75
5,3.75
0,5
1.25,5
2.5,5
3.75,5
5,5];

x=xy(:,1);
y=xy(:,2);

nebc=27;
iebc=[2;3;6;9;12;15;13;17;28;30;32;43;45;47;58;60;61;62;64;65;67;68;70;71;73;74;75];
vebc=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

nnbc=0;
inbc=[0];
vnbc=[0];

E=1*10^7;  NU=0.25; h=0.1;


fc=1;
C=E/(1-NU^2)*[1,NU,0;NU,1,0;0,0,(1-NU)/2];   
K=5/6;
A=C(3,3)*h;
D=(E*h^3)/(12*(1-NU^2)); 
