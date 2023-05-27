
clc; clear;

%Matlab code for EE768 Assignment-2 question6 

%Q3-Part1
%display('Generalized matrices for line segment')
Zeq_L = 0.75*[0.4576+1.0780i 0.1560+0.5017i 0.1535+0.3849i; 0.1560+0.5017i 0.4666+1.0482i 0.1580+0.4236i; 0.1535+0.3849i 0.1580+0.4236i 0.4615+1.0651i];
sizeof_Zeq_L = size(Zeq_L,1);
aL = eye(sizeof_Zeq_L);
bL = Zeq_L;
cL = zeros(3,3);
dL = aL;
AL = inv(aL);
BL = AL*bL;

%display('Generalized matrices for transformer')
Zbase = ((4.16/(3^(0.5)))^2)*1000/5000;
Zt_low = (0.015+0.08i)*Zbase;
dia = [Zt_low Zt_low Zt_low];
Zt_abc = diag(dia);
nt = 69*(3^(0.5))/4.16;
at = (-nt/3)*[0 2 1; 1 0 2; 2 1 0];
bt = (-nt/3)*[0 2*Zt_low Zt_low; Zt_low 0 2*Zt_low; 2*Zt_low Zt_low 0];
ct = zeros(3,3);
dt = (1/nt)*[1 -1 0; 0 1 -1; -1 0 1];
At = (1/nt)*[1 0 -1; -1 1 0; 0 -1 1];
Bt = Zt_abc;

Zc = ((4160/(3^(0.5)))^2)/300000;
%display('Generalized matrices for capacitor bank')
ac = [1 0 0; 0 1 0; 0 0 1];
bc = zeros(3,3);
dia1 = [-Zc*i -Zc*i -Zc*i]; 
Zcap = diag(dia1);
cc = inv(Zcap);
dc = ac;

ABCD =[aL bL; cL dL]*[ac bc; cc dc];

a = ABCD(1:3,1:3);
b = ABCD(1:3,4:6);
c = ABCD(4:6,1:3);
d = ABCD(4:6,4:6);
A = inv(a);
B = A*b;

%------------------------------
Zbase1 = ((0.24^2)*1000)/500;
Zt_ab = (0.009+0.03i)*Zbase1;
Zbase2 = ((0.24^2)*1000)/167;
Zt_bc = (0.01+0.016i)*Zbase2;
Zt_ca = Zt_bc;
dia2 = [Zt_ab Zt_bc Zt_ca];
Zt1_abc = diag(dia2);
nt1 = 2400/(240);
at1 = nt1*[1 -1 0; 0 1 -1; -1 0 -1];
bt1 = (nt1/3)*[Zt_ab -Zt_ab 0; Zt_bc 2*Zt_bc 0; -2*Zt_ca -Zt_ca 0];
ct1 = zeros(3,3);
dt1 = (1/(3*nt1))*[1 -1 0; 1 2 0; -2 -1 0];
At1 = (1/(3*nt1))*[2 1 0; 0 2 1; 1 0 2];
Bt1 = (1/9)*[2*Zt_ab+Zt_bc 2*Zt_bc-2*Zt_ab 0; 2*Zt_bc-Zt_ca 4*Zt_bc-Zt_ca 0; Zt_ab-4*Zt_ca -Zt_ab-2*Zt_ca 0];
%----------------------------------

% transformer taps at neutral setting and capacitor bank at load bus
rad = 3.14/180;
ELL = [(69000*cos(30*rad))+i*(69000*sin(30*rad)); (69000*cos(90*rad))-i*(69000*sin(90*rad)); (69000*cos(150*rad))+i*(69000*sin(150*rad))];
ELN_mag = 69000/(3^(1/2));
ELN = [(ELN_mag*cos(0*rad))+i*(ELN_mag*sin(0*rad)); (ELN_mag*cos(120*rad))-i*(ELN_mag*sin(120*rad)); (ELN_mag*cos(120*rad))+i*(ELN_mag*sin(120*rad))];

S3 = [(750*cos(31.79*rad))+i*(750*sin(31.79*rad)); (500*cos(25.84*rad))+i*(500*sin(25.84*rad)); (850*cos(18.19*rad))+i*(850*sin(18.19*rad))];
S4 = [(400*cos(25.84*rad))+i*(400*sin(25.84*rad)); (150*cos(36.87*rad))+i*(150*sin(36.87*rad)); (150*cos(36.87*rad))+i*(150*sin(36.87*rad))];

%Initializations
Iabc = zeros(1,3)';
IABC = zeros(1,3)';
IL = zeros(1,3)';
Vold = zeros(1,3)';
error_max = 1;

n=0;
while(error_max>0.001)
%for j = 1:50
n=n+1 
%Forward sweep
V1 = At*ELN-Bt*(Iabc+IL);
V2 = V1;
V3 = A*V2-B*(Iabc+IL);
V4 = At1*V3-Bt1*(IL);

error = abs(V4-Vold)/240;
error_max =max(error)
if error_max<=0.001
    break
end

for k = 1:3
Iabc1(k) = conj(S3(k)*1000/V3(k));
IL1(k) = conj(S4(k)*1000/V4(k));
end
Iabc = Iabc1.';
IL = IL1.';

%Backward sweep
Vold = V4;
IABC = dt*(Iabc + IL);
end
display('Line-to-ground load voltages with capacitor bank but neutral tap setting of transformer')
iteration_without_regulator = n-1
error = error_max
V3_cal = V3*120*(3^(0.5))/4160;
V3_120_unregulated = abs(V3_cal)
V4_cal = V4*120/240;
V4_120_unregulated = abs(V4_cal)


