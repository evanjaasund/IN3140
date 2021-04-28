a = [270, 30, -45]
v = [0.1; 
    -0.05; 
    -0.05]
jacobi(a, v)

function answer = jacobi(joint_angles, joint_velocities)
    new = zeros(3, 1)
    for x=1 : 3
        this = deg2rad(joint_angles(x))
        new(x) = this
    end
    t1 = new(1);
    t2 = new(2);
    t3 = new(3);
    a1 = 100.9;
    a2 = 222.1;
    a3 = 136.2;
    Jv = [-sin(t1)*((a2*cos(t2))+(a3*cos(t2+t3))) -cos(t1)*((a2*sin(t2))+(a3*sin(t2+t3))) -cos(t1)*((a3*sin(t2+t3))); 
            cos(t1)*((a2*cos(t2))+(a3*cos(t2+t3))) -sin(t1)*((a2*sin(t2))+(a3*sin(t2+t3))) -sin(t1)*((a3*sin(t2+t3)));
            0 (a2*cos(t2))+(a3*cos(t2+t3)) a3*cos(t2+t3)];
    answer = Jv*joint_velocities;
end 


