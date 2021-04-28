%Making the different symbols as we dont know the value.
syms L1 L2 L3       
syms theta_2 theta_3
g = [0; 0; -9.81];
m1 = 0.3833;
m2 = 0.2724;
m3 = 0.1406;

%Calculating the potential energy for each of the three links.
%P1
rc1 = [0; 0; L1/2];
P1 = m1*transpose(g)*rc1;

%P2
rc2 = [cos(theta_2)*L2/3; 0; ((sin(theta_2)*L2/2) + L1)];
P2 = m2*transpose(g)*rc2;

%P3
rc3 = [cos(theta_2)*L2/ + cos(theta_3)*L3/2; 0; ((sin(theta_2)*L2 + L1) + cos(theta_3)*L3/2)];
P3 = m3*transpose(g)*rc3;

%Simplifies the answers achieved from above.
P1_simplified = vpa(P1, 4)
P2_simplified = vpa(P2, 4)
P3_simplified = vpa(P3, 4)

%Computing the full potential energy.
Total_potential_energy = vpa(P1_simplified + P2_simplified + P3_simplified, 4)

