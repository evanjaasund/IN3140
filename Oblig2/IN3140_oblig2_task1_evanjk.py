import matplotlib.pyplot as plt
import numpy as np
import math
from numpy.core.numeric import ones_like
from numpy.core.shape_base import atleast_2d

#Task 1a
def forward(joint_angles): #function that takes joins angles and return cartesian coordinates

    c1 = np.cos(np.deg2rad(joint_angles[0]))
    c2 = np.cos(np.deg2rad(joint_angles[1]))
    c3 = np.cos(np.deg2rad(joint_angles[2]))
    s1 = np.sin(np.deg2rad(joint_angles[0]))
    s2 = np.sin(np.deg2rad(joint_angles[1]))
    s3 = np.sin(np.deg2rad(joint_angles[2]))
    d1 = 100.9; a2 = 222.1; a3 = 136.2

    a_0_3 = np.array([[(a3*c1*c2*c3)-(a2*c1*s2*s3)+(a2*c1*c2)], [(a3*s1*c2*c3)-(a3*s1*s2*s3)+(a2*s1*c2)], [(a3*s2*c3)+(a3*c2*s3)+(a2*s2)+(d1)], [1]])

    # to_task_coord = np.array([[0, 1, 0, 750], [-1, 0, 0, 250], [0, 0, 1, -100], [0, 0, 0, 1]])
    # a = to_task_coord@a_0_3

    x_value = np.round(a_0_3[0], 4)
    y_value = np.round(a_0_3[1], 4)
    z_value = np.round(a_0_3[2], 4)

    cart_coord = [x_value[0], y_value[0], z_value[0]]
    return cart_coord

#Task1b & 1d
def inverse(cart_coord): #function that takes cartesian coordinates as an argument and returns the joint angles

    d1 = 100.9; a2 = 222.1; a3 = 136.2
    solutions = []

    r1 = math.sqrt((cart_coord[0]**2)+(cart_coord[1]**2))
    r2 = cart_coord[2]-d1
    r3 = math.sqrt((r1**2)+(r2**2))
    my1 = math.acos((r3**2-a2**2-a3**2)/-(2*a2*a3))
    my2 = math.acos((a3**2-r3**2-a2**2)/-(2*r3*a2))
    my3 = math.atan2(r2, r1)

    #elbow down & Task1c
    t1 = round(np.rad2deg(math.atan2(cart_coord[1], cart_coord[0])), 4)
    t2 = round(np.rad2deg((my3 - my2)), 4)
    t3 = round(np.rad2deg((math.pi - my1)), 4)
    solutions.append([t1, t2, t3])

    #elbow down rotated
    t4 = -t1
    t5 = 180 - round(np.rad2deg((my3 - my2)), 4)
    t6 = t3
    solutions.append([t4, t5, t6])

    #elbow up 
    t7 = t1
    t8 = round(np.rad2deg((my3 + my2)), 4)
    t9 = -t3
    solutions.append([t7, t8, t9])

    #elbow up rotated 
    t10 = -t1
    t11 = 180 - round(np.rad2deg((my3 + my2)), 4)
    t12 = -t3
    solutions.append([t10, t11, t12])

    return solutions


joint_angles = [270, 30, -45]  
print("Joint angles sent in: ")
cartesian_coordinates = forward(joint_angles) #cart_coord assigned the joint angles
print(joint_angles)                           #sending the joint angles in
print(" ")

print("Corresponding coordinates: ")
print(cartesian_coordinates)
solutions = inverse(cartesian_coordinates)    #solutions assigned solution array
print(" ")


print("Joint angles returned from inverse function: \n" )
print(solutions[2])                           #returns -90 for theta1, which is the same as 270

