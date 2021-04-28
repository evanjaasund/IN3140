%%Task 2b)
syms I_1x I_1y I_1z I_2x I_2y I_2z I_3x I_3y I_3z
syms L1 L2 L3       
syms theta_1 theta_2 theta_3
syms q_dotx q_doty q_dotz q_ddotx q_ddoty q_ddotz
g = [0; 0; -9.81];
thetas = [theta_1, theta_2, theta_3];
q_dot = [q_dotx; q_doty; q_dotz];
q_ddot = [q_ddotx; q_ddoty; q_ddotz];

m1 = 0.3833;
m2 = 0.2724;
m3 = 0.1406;

%Calculating the potential energy for each of the three links.
%P1
rc1 = [0; 0; L1/2]
P1 = m1*transpose(g)*rc1;

%P2
rc2 = [cos(theta_2)*L2/3; 0; ((sin(theta_2)*L2/2) + L1)]
P2 = m2*transpose(g)*rc2;

%P3
rc3 = [cos(theta_2)*L2/ + cos(theta_3)*L3/2; 0; ((sin(theta_2)*L2 + L1) + cos(theta_3)*L3/2)]
P3 = m3*transpose(g)*rc3;

%Simplifies the answers achieved from above.
P1_simplified = vpa(P1, 4);
P2_simplified = vpa(P2, 4);
P3_simplified = vpa(P3, 4);

%Computing the full potential energy.
Total_potential_energy = vpa(P1_simplified + P2_simplified + P3_simplified, 4);

I1 = [I_1x 0 0 
    0 I_1y 0 
    0 0 I_1z];

I2 = [I_2x 0 0
    0 I_2y 0
    0 0 I_2z];

I3 = [I_3x 0 0
    0 I_3y 0
    0 0 I_3z];


%Rotation matrices.
R1 = [1 0 0;
    0 1 0;
    0 0 1];
R2 = [cos(theta_2) -sin(theta_2) 0
    sin(theta_2) cos(theta_2) 0
    0 0 1];
R3 = [cos(theta_3) -sin(theta_3) 0
    sin(theta_3) cos(theta_3) 0
    0 0 1];

%Making the tranformation matrices.
R01 = R1;
R02 = R01 * R2;
R03 = R02 * R3;

%Filling in the jacobian.      
Jv1 = [0 0 0; 0 0 0; 0 0 0];
Jv2 = [-sin(theta_1)*(L2*cos(theta_2))  -cos(theta_1)*(L2*sin(theta_2)) 0; -cos(theta_1)*(L2*cos(theta_2)) -sin(theta_1)*(L2*sin(theta_2)) 0; 0 L2*cos(theta_2) 0];
Jv3 = [-sin(theta_1)*(L2*cos(theta_2)) + L3*cos(theta_2)*cos(theta_3) -cos(theta_1)*(L2*sin(theta_2)) + L3*sin(theta_2)*sin(theta_3) -cos(theta_1)*(L3*sin(theta_2)*sin(theta_3));
    -cos(theta_1)*(L2*cos(theta_2)) + (L3*cos(theta_2)*cos(theta_3)) -sin(theta_1)*(L2*sin(theta_2)) + L3*sin(theta_2)*sin(theta_3) -sin(theta_1)*(L3*sin(theta_2))*sin(theta_3);
    0 L2*cos(theta_2) + L3*cos(theta_2)*cos(theta_3) L3*cos(theta_2)*cos(theta_3)];

Jw1 = [0 0 0; 0 0 0; 1 0 0];
Jw2 = [0 sin(theta_1) 0; 0 -cos(theta_1) 0; 1 0 0];
Jw3 = [0 sin(theta_1) sin(theta_1);
    0 -cos(theta_1) -cos(theta_1);
    1 0 0];


Jv1_T = transpose(Jv1);
Jw1_T = transpose(Jw1);
K1 = m1*Jv1_T* Jv1 + Jw1_T * R01 * I1 * transpose(R01) * Jw1;

Jv2_T = transpose(Jv2);
Jw2_T = transpose(Jw2);
K2 = m2*Jv2_T* Jv2 + Jw2_T * R02 * I2 * transpose(R02) * Jw2;

Jv3_T = transpose(Jv3);
Jw3_T = transpose(Jw3);
K3 = m3*Jv3_T* Jv3 + Jw3_T * R03 * I3 * transpose(R03) * Jw3;

K = 1/2 * transpose(q_dot) * (K1 + K2 + K3) * q_dot;

%% Task2d-e

%From equation 14 and 15 in the mandatory assignment we can see that 
%the K made above can is the same as D(q). We therefore already have the
%first term. Therefore we can say that D:
D = K1 + K2 + K3;

%How to compute g_k(q) is to be found in the 41st slide of the dynamics-PDF.
%We then need to partial derivate the potential energy for each link with 
%the perspective of theta as we have no q.


g1 = diff(P1, theta_1);
g2 = diff(P2, theta_2);
g3 = diff(P3, theta_3);

G = [g1; g2; g3];

C = getC(D, thetas);

T = D * q_ddot + C * q_dot + G


%To find the last part of the equation (C) we use the formula found at page
%41 in the dynamics-PDF.
function matrixC = getC(D, syms)
    matrixC = sym(zeros(3,3));
    for k=1:3
        for j=1:3
            for i=1:3
                term1 = diff(D(k,j), syms(i));
                term2 = diff(D(k,i), syms(j));
                term3 = diff(D(i,j), syms(k));
                matrixC(k,j) = matrixC(k,j) + 1/2 * (term1+term2-term3);
            end
        end
    end
end



