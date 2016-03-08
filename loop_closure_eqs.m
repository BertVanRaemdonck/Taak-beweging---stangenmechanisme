%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F=loop_closure_eqs(phi_init,phi2,r2l,r2k,r3,a,b,r6l,r6k,r7,r8l,r8k,r10,r11,r12,x4,y4,x7,y7,y9)

% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% declaration of new variable
r6 = r6k + r6l;

% copy initial values of unknown angles phi3 and phi4
phi3=phi_init(1);
phi4=phi_init(2);
x5=phi_init(3);
phi6=phi_init(4);
phi7=phi_init(5);
phi8=phi_init(6);
x9=phi_init(7);
phi10=phi_init(8);
x11=phi_init(9);
phi12=phi_init(10);

% loop closure equations:
F(1)=-r2l*cos(phi2)+r12*cos(phi12)-(r2l+r12-x11);
F(2)=-r2l*sin(phi2)+r12*sin(phi12);
F(3)=-r10*cos(phi10)+r8l*cos(phi8)-x9-x11;
F(4)=-r11-r10*sin(phi10)+r8l*sin(phi8)-y9;
F(5)=-r7*cos(phi7)+r6k*cos(phi6)+x5*cos(phi4)-x4+x7;
F(6)=-r7*sin(phi7)+r6k*sin(phi6)+x5*sin(phi4)-y4+y7;
F(7)=-r7*cos(phi7)+r6*cos(phi6)-r8k*cos(phi8)-(r2l+r12+x9)+x7;
F(8)=-r7*sin(phi7)+r6*sin(phi6)-r8k*sin(phi8)-y9+y7;
F(9)=r2k*cos(phi2-pi/2)-r3*cos(phi3)+a*cos(phi4)-b*cos(phi4+pi/2)-x4;
F(10)=r2k*sin(phi2-pi/2)-r3*sin(phi3)+a*sin(phi4)-b*sin(phi4+pi/2)-y4;

