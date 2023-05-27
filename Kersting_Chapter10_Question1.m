clc; clear;

%Matlab code for EE768 Assignment-2 question-3 
%Developed by Divya M 18V972020

%Part1
display('Generalized matrices for line segment')
Zeq_L = 0.75*[0.4576+1.0780i 0.1560+0.5017i 0.1535+0.3849i; 0.1560+0.5017i 0.4666+1.0482i 0.1580+0.4236i; 0.1535+0.3849i 0.1580+0.4236i 0.4615+1.0651i];
sizeof_Zeq_L = size(Zeq_L,1);
aL = eye(sizeof_Zeq_L)
bL = Zeq_L
cL = zeros(3,3)
dL = aL
AL = inv(aL)
BL = AL*bL

display('Generalized matrices for transformer')
Zbase = ((4.16^2)*1000)/5000;
Zt_low = (0.015+0.08i)*Zbase;
dia = [Zt_low Zt_low Zt_low];
Zt_abc = diag(dia);
nt = (69*sqrt(3))/4.16;
at = (-nt/3)*[0 2 1; 1 0 2; 2 1 0]
bt = (-nt/3)*[0 2*Zt_low Zt_low; Zt_low 0 2*Zt_low; 2*Zt_low Zt_low 0]
ct = zeros(3,3)
dt = (1/nt)*[1 -1 0; 0 1 -1; -1 0 1]
At = (1/nt)*[1 0 -1; -1 1 0; 0 -1 1]
Bt = Zt_abc

%Part2
rad = 3.14/180;
ELL = [(69000*cos(30*rad))+i*(69000*sin(30*rad)); (69000*cos(90*rad))-i*(69000*sin(90*rad)); (69000*cos(150*rad))+i*(69000*sin(150*rad))];
ELN_mag = 69000/(sqrt(3));
ELN = [(ELN_mag*cos(0*rad))+i*(ELN_mag*sin(0*rad)); (ELN_mag*cos(120*rad))-i*(ELN_mag*sin(120*rad)); (ELN_mag*cos(120*rad))+i*(ELN_mag*sin(120*rad))];

S3 = [(750*cos(31.79*rad))+i*(750*sin(31.79*rad)); (500*cos(25.84*rad))+i*(500*sin(25.84*rad)); (850*cos(18.19*rad))+i*(850*sin(18.19*rad))];

%Initializations
Iabc = zeros(1,3)'; %current at node3
IABC = zeros(1,3)'; %source current
Vold = zeros(1,3)';
error_max = 1;

n=0;
while(error_max>0.001)  %tolerance is 0.001
  n=n+1 ; 
%Forward sweep
V1 = At*ELN-Bt*Iabc;    
V2 = V1;
V3 = AL*V2-BL*Iabc;

error = abs(V3-Vold)*sqrt(3)/4160;
error_max =max(error);
if error_max<=0.001
    break
end

for i = 1:3
Iabc1(i) = conj(S3(i)*1000/V3(i));
end
Iabc = Iabc1.';

%Backward sweep
Vold = V3;
IABC = dt*Iabc;
end

display('Results of system without voltage regulators and shunt capacitors')
total_iterations = n-1
error = error_max
[th,mag] =cart2pol(real(V3),imag(V3));
th_deg = th*180/3.14;
V3_cal = V3*120*sqrt(3)/4160;
V3_120 = abs(V3_cal)
V3_phase =mag