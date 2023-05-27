
clc; clear;

%Matlab code for EE768 Assignment-2 question5 
%Developed by Divya M

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
Zbase = ((4.16^2)*1000)/5000;
Zt_low = (0.015+0.08i)*Zbase;
dia = [Zt_low Zt_low Zt_low];
Zt_abc = diag(dia);
nt = 69*sqrt(3)/4.16;
at = (-nt/3)*[0 2 1; 1 0 2; 2 1 0];
bt = (-nt/3)*[0 2*Zt_low Zt_low; Zt_low 0 2*Zt_low; 2*Zt_low Zt_low 0];
ct = zeros(3,3);
dt = (1/nt)*[1 -1 0; 0 1 -1; -1 0 1];
At = (1/nt)*[1 0 -1; -1 1 0; 0 -1 1];
Bt = Zt_abc;

Zc = ((4160/sqrt(3))^2)/300000;
%display('Generalized matrices for capacitor bank')
ac = [1 0 0; 0 1 0; 0 0 1];
bc = zeros(3,3);
dia1 = [1*i/Zc;1*i/Zc; 1*i/Zc ]; 
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

% transformer taps at neutral setting and capacitor bank at load bus
rad = 3.14/180;
ELL = [(69000*cos(30*rad))+i*(69000*sin(30*rad)); (69000*cos(90*rad))-i*(69000*sin(90*rad)); (69000*cos(150*rad))+i*(69000*sin(150*rad))];
ELN_mag = 69000/(sqrt(3));
ELN = [(ELN_mag*cos(0*rad))+i*(ELN_mag*sin(0*rad)); (ELN_mag*cos(120*rad))-i*(ELN_mag*sin(120*rad)); (ELN_mag*cos(120*rad))+i*(ELN_mag*sin(120*rad))];

S3 = [(750*cos(31.79*rad))+i*(750*sin(31.79*rad)); (500*cos(25.84*rad))+i*(500*sin(25.84*rad)); (850*cos(18.19*rad))+i*(850*sin(18.19*rad))];

%Initializations
Iabc = zeros(1,3)';
IABC = zeros(1,3)';
Vold = zeros(1,3)';
error_max = 1;

n=0;
while(error_max>0.001)
n=n+1 ; 
%Forward sweep
V1 = At*ELN-Bt*(Iabc);
V2 = V1;
V3 = A*V2-B*(Iabc);

error = abs(V3-Vold)*sqrt(3)/4160;
error_max =max(error);
if error_max<=0.001
    break
end

for k = 1:3
Iabc1(k) = conj(S3(k)*1000/V3(k));

end
Iabc = Iabc1.';

%Backward sweep
Vold = V3;
IABC = dt*(Iabc);
end
display('Line-to-ground load voltages with capacitor bank but neutral tap setting of transformer')
iteration_without_regulator = n-1
error = error_max
V3_cal = V3*120/4160;
V3_120_unregulated = abs(V3_cal)

% Transformer tap settings
PT = 20;
CTp = 500;
CTs = 5;
CT = CTp/CTs;
Tap = [0; 0; 0];

n1=0;
ans = 0;
while(ans==0)
 n1=n1+1;   
for k =1:3
Zeq1(k) = (V2(k)-V3(k))/(Iabc(k));
end
Zeq = Zeq1.';
Zavg = (1/3)*sum(Zeq);
Z_volt = Zavg*CTp/PT;
Z_ohm = Z_volt/CTs;
V_reg = V2/((3^0.5)*PT);
I_com = Iabc/((3^0.5)*CT);
V_relay = V_reg-Z_ohm*I_com;

dia = [1-0.00625*Tap(1) 1-0.00625*Tap(2) 1-0.00625*Tap(3)];

%display('Generalized matrices for voltage regulator')
ar = diag(dia);
br =zeros(3,3);
cr = br;
dr = inv(ar);
Ar = dr;
Br = br;

%Initializations
Iabc = zeros(1,3)';
IABC = zeros(1,3)';
Ir = zeros(1,3)';   
Vold = zeros(1,3)';
error_max = 1;

while(error_max>0.001)
  
%Forward sweep
V1 = At*ELN-Bt*Ir;
V2 = Ar*V1-Br*(Iabc);
V3 = A*V2-B*(Iabc);

error = abs(V3-Vold)/2400;
error_max =max(error);
if error_max<=0.001
    break
end

for k = 1:3
Iabc1(k) = conj(S3(k)*1000/V3(k));
end
Iabc = Iabc1.';


%Backward sweep
Vold = V3;
Ir = dr*(Iabc);       
IABC = dt*Ir ;     
end
V3_120 = V3*120/4160;
n = [1; 1; 1];
for k =1:3
   if (abs(V3_120(k))<121)
        Tap(k) = Tap(k)+1;
        %encounter
        n(k) = 0;
    end
end
%Tap
ans = min(n);
end

display('Required tap settings')
Tap
display('Relay and load voltages under regulated conditions')
iteration = n1-1
V3_regulated = abs(V3_120)

