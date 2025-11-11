
%Project: Fixed-Wing Hybrid-Electric //can be altered.

clc
clear
close all

% BODY NX=7.0,ITYPE=1.0,METHOD=1.0, [ft]
XB=[0.00, 5.95, 20.44, 34.93, 43.79, 50.38, 52.00];
% XB= [];
% for ii=1:7
%     XB(ii)= (ii-1)*(36/7);
% end
ZU=[0.00, 2.67, 6.70, 5.64, 4.11, 3.72, 3.62] ;
ZL=[0.00, 1.73, 2.54, 1.05, -1.07, -2.91, 3.52];
RB=[0.00, 4.78, 8.34, 6.38, 3.27, 0.77, 0.00];

% Mach= 0.287; 
Weight= 15000;              %[lb]
ALT=   15000.0;             %[ft]
Velocity= 400.0;            %[knots]
Velocity= Velocity*1.69;    %[ft/s]
Pressure= 1.195*10^3;       %[lb/ft^2]
Temp= 465.2;                %[deg R]
R= 1716.46;                 %[ft-lb/slug-deg R]
Density= Pressure/(R*Temp); %[lb/ft^3]
xcg_tilda_cg= 0.25;         %cg location
v= 0.2293*10^(-3);          %[ft^2/s]   from A.2
Q= 0.5*Density*Velocity^2;
sos= sqrt(1.4*R*Temp);      
Mach= Velocity/sos;
beta= sqrt(1-Mach^2);       %Mach parameter

% SYNTHS parameters
XW.SYN=     16.78;      %nose to wing LE
XW.SYN=     16.78;
ZW.SYN=    -1.97;       
ALIW.SYN=   3.00;       %wing incedence angle
XH.SYN=    42.99;       %nose to HT LE
ZH.SYN=     2.74;
ALIH.SYN=   0.00;
XV.SYN=    42.99;       %nose to VT LE (equal to HT for V-tail)
ZV.SYN=     2.74; 
XP=         9.00;       %nose to props


% WING parameters
CHRDR.w=   8.75;        %c_r
CHRDTP.w=  5.0;         %c_t
SSPN.w=   60.0/2;       %b_w/2 semi-span
SSPNE.w=  60.0/2;       %b_w*/2 exposed semi-span
SAVSI.w=   1.0;         %sweep LE
CHSTAT.w=  0.0;         % of mac sweep angle reference
TWISTA.w=  0.0;         %twist angle
DHDADI.w=  3.0;         %dihedral
% TC.w=      0.18;


% HORIZONTAL TAIL parameters
CHRDR.ht=   5.94;
CHRDTP.ht=  3.405;
SSPN.ht=    (24.0/2)*sind(45);
SSPN.ht=    25.279;
SSPNE.ht=   (24.0/2-4.1)*sind(45);
SSPNE.ht=   25.279;
SAVSI.ht=  25.0;
CHSTAT.ht = 0.0;    
TWISTA.ht = 0.0;    
DHDADI.ht=  0.0;    
% TC.ht =     0.12;


% VERTICAL TAIL parameters
CHRDR.vt=   CHRDR.ht;
CHRDTP.vt=  CHRDTP.ht;
SSPN.vt=    SSPN.ht;
SAVSI.vt=   SAVSI.ht;
SSPNE.vt =  SSPNE.ht;
CHSTAT.vt = 0;


% %VERTICAL TAIL RH OUTTER TOP
% XV.vt_2t=      29.65;
% YV.vt_2t=       7.487;
% ZV.vt_2t=       2.3;
% CHRDR.vt_2t=    4.023;
% CHRDBP.vt_2t=   0.0;
% CHRDTP.vt_2t=   2.30;
% SSPN.vt_2t=     3.595;
% SSPNOP.vt_2t=   0.0;
% SAVSO.vt_2t=    0.0;
% SAVSI.vt_2t=    0.0;
% CHSTAT.vt_2t=   0.58;
% VERTUP.vt_2t=      1;
% % TC.vt=          0.12;
% 
% 
% %VERTICAL TAIL LH OUTTER TOP
% XV.vt_3t=         29.65;
% YV.vt_3t=         -7.487;
% ZV.vt_3t=          2.3;
% CHRDR.vt_3t=    4.023;
% CHRDBP.vt_3t=   0.0;
% CHRDTP.vt_3t=   2.30;
% SSPN.vt_3t=     3.595;
% SSPNOP.vt_3t=   0.0;
% SAVSO.vt_3t=    0.0;
% SAVSI.vt_3t=    0.0;
% CHSTAT.vt_3t=   0.58;
% VERTUP.vt_3t=      1;
% % TC.vt=          0.12;
% 
% 
% %VERTICAL TAIL LH OUTTER BOTTOM
% XV.vt_3b =        29.65;
% YV.vt_3b =        -7.487;
% ZV.vt_3b =         2.3;
% CHRDR.vt_3b =   4.023;
% CHRDBP.vt_3b =  0.0;
% CHRDTP.vt_3b =  2.30;
% SSPN.vt_3b =    1.822;
% SSPNOP.vt_3b =  0.0;
% SAVSO.vt_3b =   0.0;
% SAVSI.vt_3b =   0.0;
% CHSTAT.vt_3b =  0.58;
% VERTUP.vt_3b =     0;
% TC.vt =         0.12;


%% Calculations

%Wing
c_r.w= CHRDR.w;     %chord root [ft]
c_t.w= CHRDTP.w;    %chord tip [ft]
Lambda.w= SAVSI.w;  %LE sweep [deg]
Gamma.w= DHDADI.w;  %dihedral [deg]
b.w= SSPN.w*2;      %span [ft]
i.w= ALIW.SYN;      %incidence angle of the wing
a_0.w= 0.106; %a_0.w= 0.110;%0009     %wing mach zero lift-curve slope [per deg] AF#2415
a_0.wrad= a_0.w*180/pi; 
k= a_0.wrad/(2*pi);    
ac.w= 0.243;        %x location of AC on wing
Cm.ac= -0.045; %Cm.ac= 0;%0009      %2D Cm_ac of wing AF#2415
A.w= b.w*2/(c_r.w+c_t.w);   %Aspect Ratio
lambda.w= c_t.w/c_r.w;      %Taper Ratio
Lambda.w_c2= atand(tand(Lambda.w)-(4*(1/2)*(1-lambda.w))/(A.w*(1+lambda.w))); %wing sweep at half-chord
a.wrad= 2*pi*A.w/(2+sqrt((beta*A.w/k)^2*(1+((tand(Lambda.w_c2)^2)/(beta^2)))+4));       %wing lift-curve slope [per rad]
a.w= a.wrad*pi/180;     %wing lift-curve slop [per deg]
c_bar.w= (2/3)*c_r.w*((1+lambda.w+lambda.w^2)/(1+lambda.w)); %mac [ft]
y_bar.w= (b.w/6)*(1+2*lambda.w)/(1+lambda.w);                %y location of mac [ft]
x_bar.w= y_bar.w*tand(Lambda.w);                             %x location of mac [ft]
S.w= b.w^2/A.w;         %wing area [ft^2]
e.w= 0.985;            %from fig. 4.22 on pg. 4.5

xcg= xcg_tilda_cg*c_bar.w;      %x-location of cg wrt LE
XCG= XW.SYN+xcg+x_bar.w;        %x-location of cg wrt nose

%Horizontal Tail
c_r.ht= CHRDR.ht;       %chord root [ft]
c_t.ht= CHRDTP.ht;      %chord tip [ft]
Lambda.ht= SAVSI.ht;    %LE sweep [deg]
b.ht= SSPN.ht*2;        %span [ft]
i.ht= ALIH.SYN;         %incidence angle of the tail
a_0.ht= 0.09596;        %HT mach zero lift-curve slope[per deg]
a.ht= 0.09949;          %wing lift-curve slope [per deg]
ac.ht= 0.25;            %x location of AC on HT
n.ht= 0.95;              %estimated
A.ht= b.ht*2/(c_r.ht+c_t.ht);   %Aspect Ratio
S.ht= b.ht^2/A.ht;              %wing area [ft^2]
lambda.ht= c_t.ht/c_r.ht;       %Taper Ratio
c_bar.ht= (2/3)*c_r.ht*((1+lambda.ht+lambda.ht^2)/(1+lambda.ht)); %mac [ft]
y_bar.ht= (b.ht/6)*(1+2*lambda.ht)/(1+lambda.ht);                %y location of mac [ft]
x_bar.ht= y_bar.ht*tand(Lambda.ht);                             %x location of mac [ft]
l.ht= (XH.SYN+x_bar.ht+ac.ht*c_bar.ht)-XCG;

%Fuselage
dZi= ZU-ZL;                         %width array of fuselage [ft]
l.f= XW.SYN+0.25*c_r.w;             %nose to wing root chord [ft]
L.f= 33.6;                          %nose to tail cone [ft]
w.f= max(dZi);                        %max width/diameter [ft]
K.f= 0.002322*exp(5.037*(l.f/L.f)); %correction factor [/deg]


%Vertical Tail
c_r.vt= CHRDR.vt;       %chord root [ft]
c_t.vt= CHRDTP.vt;      %chord tip [ft]
b.vt= SSPN.vt+SSPN.vt;    %total span [ft]
b.vt_prime= SSPN.vt;         %upper span [ft]
lambda.vt= c_t.vt/c_r.vt;       %Taper Ratio
c_bar.vt= (2/3)*c_r.vt*((1+lambda.vt+lambda.vt^2)/(1+lambda.vt)); %mac [ft]
y_bar.vt= (b.vt/6)*(1+2*lambda.vt)/(1+lambda.vt);                %y location of mac [ft]
x_bar.vt= y_bar.vt*tand(Lambda.ht);                             %x location of mac [ft]
bvp_bv= b.vt_prime/b.vt;        %ratio to find Aeff
Aeff_Avt= 1.55;                 %ratio from Fig. 5.3.1.1-24a on pg. 9.13
A.vt= b.vt*2/(c_r.vt+c_t.vt);   %Aspect Ratio
A.eff= Aeff_Avt*A.vt;           %effective aspect ratio
t_c= 12;                        
angle_TE= 20;                   %approx., use to find CYB.veff
CYB.veff= 3.25;                 %CYB from Fig 5.3.1.1-24b on pg 9.16 [/rad]
CYB.veffrad= 2.5;               %alternate
CYB.veff= CYB.veffrad*pi/180;   %alternate in deg
S.vt= b.vt^2/A.vt;              %wing area [ft^2]
ratio.c= w.f/b.vt;              %inital ratio not sure what its called
intercept.left= b.ht/L.f;       %b_HT / length fuselage
CYB.vWBH_veff= 0.68;            %ratio from Fig 5.3.1.1-24c on pg 9.16
CYB.vWBH_veff= 0.725;
CYB.tvt= -(CYB.veff/57.3)*CYB.vWBH_veff*2*S.vt/S.w;

alpha.ZL= -2.0; %alpha.ZL= 0;%0009         %alpha zero-lift from Table G.1 on pg. A.27
                   

%Side force gradient
depsilon_dalpha= 2*(a.wrad)/(pi*e.w*A.w); %theoretical downwash gradient

%Calculating CLs
a_0_M= a_0.w/beta;  % including mach effects
CL.o= a_0_M*(i.w-alpha.ZL);   %C_L_0
CL.crit= Weight/(0.5*Density*Velocity^2*S.w);   
CL.alpha= a.w+a.ht*n.ht*(S.ht/S.w)*(1-depsilon_dalpha);
CL.cruisewinc= CL.alpha * i.w; % Additional lift from incidence angle
alpha.crit= (CL.crit-(CL.o-a_0_M*i.w))/CL.alpha; %[deg]
CL.crit= CL.crit+CL.cruisewinc;


%Calculating CMs
CM.alpha_f= K.f*w.f^2*L.f/(S.w*c_bar.w); %[/deg]
CM.alpha= a.w*(xcg_tilda_cg-ac.w)+CM.alpha_f-a.ht*n.ht*(S.ht*l.ht/(S.w*c_bar.w))*(1-depsilon_dalpha);
% CM.o= -CM.alpha*alpha.crit; %NOT CORRECT
CM.of= 0;           %assume no fuselage moment
CM.acw= -0.045; %CM.acw=0;%0009    %Table G.1 on pg A.27
CM.o= CM.of + CM.acw + a.w*(i.w-alpha.ZL)*(xcg_tilda_cg-ac.w) - a.ht*n.ht*(S.ht*l.ht/(S.w*c_bar.w))*(i.ht-depsilon_dalpha*(i.w-alpha.ZL));

%Remaining CYB Calcs
%Wing-Body:
z= 0.95;  %[ft]
d_atz= 5.5; %[ft]
z_d2= z/(d_atz/2);
K.wb= 1.175;        %correction factor from Fig 5.2.1.1-7 on pg 9.6
fineness_ratio= L.f/w.f; %ratio for k2-k1
K2_K1= 0.875;       %mass factor from Fig 4.2.1.1-20a on pg 9.7
x.l= 14.0;        %station loc at start of decrease
x.o= (0.378+0.527*(x.l/l.f))*l.f; %station for c-sec
S.o= pi*2.45^2;     %cross-sectional area where flow ceases to be potential
CYB.wb= -(2/57.3)*K.wb*K2_K1*(S.o/S.w)-0.0001*abs(Gamma.w); 

%Props
bhp= 450;
n.pe= 2;        %number of props
D.p= 8.25;      %model info
T.c= 550*bhp*n.pe/(Density*Velocity^3*D.p^2);
f=1+0.82*log(1+T.c);
beta_pitch= (29-14)*0.75+14;
nom_pitch= 8.28;
S.p= pi*(nom_pitch/2)^2;
CYB.Tc0= 0.120;     %fig. 4.6.1-25a [per rad]
CYB.Tc0= CYB.Tc0*pi/180;     %fig. 4.6.1-25a [per deg]
CYB.prop= -(n.pe/57.3)*f*CYB.Tc0*(S.p/S.w);

CYB.total= CYB.wb + CYB.tvt +CYB.prop;

%Calculating CNBs
CNB.wd= -CL.crit*(0.075/(57.3)^2)*Gamma.w; %wing due to dihedral

Lambda.w_c4= atand(tand(Lambda.w)-(4*(1/4)*(1-lambda.w))/(A.w*(1+lambda.w))); %wing sweep at quater-chord

CNB.wsw_M0= (CL.crit^2/57.3)*((1/(4*pi*A.w))- ...
    (tand(Lambda.w_c4)/(pi*A.w*(A.w+4*cosd(Lambda.w_c4))))* ...
    (cosd(Lambda.w_c4)-A.w/2-A.w^2/(8*cosd(Lambda.w_c4))-(6/A.w)* ...
    (xcg_tilda_cg-ac.w)*sind(Lambda.w_c4))); %wing due to sweep at M=0

B = sqrt(1-(Mach^2*cosd(Lambda.w_c4)^2));

CNB.wsw_term1 = (A.w+4*cosd(Lambda.w_c4))/(A.w*B+4*cosd(Lambda.w_c4)); 

CNB.wsw_term2= (A.w^2*B^2+4*A.w*B*cosd(Lambda.w_c4)-8*cosd(Lambda.w_c4)^2)/(A.w^2+4*A.w*cosd(Lambda.w_c4)-8*cosd(Lambda.w_c4)^2);

CNB.wsw= CNB.wsw_term1 * CNB.wsw_term2 *CNB.wsw_M0; %wing due to sweep

%Calculate ratios to find Kn
fig8.ratio1= XCG/L.f;
S.bs= 0; %initial body side area variable
nx= length(XB);
for int=1:(nx-1)
    j= (1/2)*(dZi(int)+dZi(int+1))*(XB(int+1)-XB(int));
    S.bs= S.bs + j; %body side area counter
end
fig8.ratio2= L.f^2/S.bs;
XB_interp=[L.f/4, 3*L.f/4]; 
h_vbl= interp1(XB, dZi, XB_interp); %width at 1/4 and 3/4 L.f
fig8.ratio3= sqrt(h_vbl(1)/h_vbl(2));
fig8.ratio4= w.f/w.f;
K.N= 0.001; %Fig. 5.2.3.1-8 on pg. 9.31 using ratios 1-4 [per deg]

Re.lf= Velocity*L.f/v;
K.rl= 0.205*log(Re.lf)-1.833; %also check Fig 5.2.3.1-9 on pg. 9.32
 
CNB.bw= -K.N*K.rl*(S.bs/S.w)*(L.f/b.w);

l.v= (XV.SYN + x_bar.vt + ac.ht * c_bar.vt ) - XCG;

CNB.tvt= -CYB.tvt*l.v/b.w;

x.p= XCG- XP;

CNB.prop= CYB.prop*x.p/b.w;
 
CNB.total= CNB.wd + CNB.wsw + CNB.bw + CNB.tvt + CNB.prop;

%Calculating ClBs
ClB_CL.sw_c2= -0.00005;     %from Fig. 5.1.2.1-27 on pg. 9.37
ClB.wsw_ratio1= Mach*cosd(Lambda.w_c2);
ClB.wsw_ratio2= A.w/cosd(Lambda.w_c2);
K.Msw= 1.025;               %from Fig. 5.1.2.1-28a
l.f_c2= XW.SYN+(1/2)*c_r.w+(b.w/2)*tand(Lambda.w_c2); %from nose to c_t/2
ClB.wsw_ratio3= l.f_c2/b.w;
K.f_c2= 1.0;                %from Fig. 5.2.2.1-26 on pg. 9.41

ClB.wsw= CL.crit*ClB_CL.sw_c2*K.Msw*K.f_c2;

ClB_CL.A= -0.0006667;      %from Fig. 5.1.2.1-28b on pg. 9.38

ClB.wA= CL.crit*ClB_CL.A;

dClB_tw_c4= -0.000033;   %Fig. 5.1.2.1-30b [per deg^2]

ClB.wtw= TWISTA.w*tand(Lambda.w_c4)*dClB_tw_c4;

ClB_Di_Di= -0.00022;    %Fig. 5.1.2.1-29b on pg. 9.39
K.Mdi= 1.15;            %Fig. 5.1.2.1-30a on pg. 9.40
d.eq= 2*interp1(XB, RB, (XW.SYN+x_bar.w));  %diameter at wing root chord

ClB.wdi= Gamma.w*ClB_Di_Di*K.Mdi - Gamma.w*(0.0005*sqrt(A.w)*(d.eq/b.w)^2);

z_prm.w= -interp1(XB, ZL, (XW.SYN+c_r.w/4));

ClB.bw= (1.2*sqrt(A.w)/57.3)*(2*d.eq*z_prm.w/(b.w)^2);

h.vt= -ZV.SYN*cosd(alpha.crit)-(XV.SYN+ac.ht-XCG)*sind(alpha.crit);

ClB.tvt= CYB.tvt*h.vt/b.w;

ClB.total= ClB.wsw + ClB.wA + ClB.wtw + ClB.wdi +ClB.bw + ClB.tvt;



%Print in window
fprintf('At Xcg= 32%% chord:\n')
%Coeff of Lifts
fprintf('\nCL_0 = %.4f\nCL_cr = %.4f\nCL_a = %.4f\n', CL.o, CL.crit, CL.alpha) 
%Coeff of Moments
fprintf('\nAoA = %.4f%c\nCM_0 = %.4f\nCM_a = %.4f\n', alpha.crit, char(176), CM.o, CM.alpha)
% %Critical AoA
% fprintf('\nAoA = %.4f%c\n', alpha.crit, char(176))
%Coeff of dynamic stabilities
fprintf('\nCYB = %.6f\nCNB = %.6f\nClB = %.6f\n', CYB.total, CNB.total, ClB.total) 
 
