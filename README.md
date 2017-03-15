# Two-D-Motion-Simulator
The classes in this repository can be used to simulate two-dimensional motion in two different formats. Both used Open Source 
Physics packages for the visual components of the simulaion. 

The first project is to simulate the motion of a ball that is subject to both gravity and air resistance. Specifically, it is
initialized to simulate the conditions of the Green Monster wall in Fenway park for the purpose of finding the optimal angle 
and smallest velocity necessary to hit a homerun over the wall. First, the program initializes an array of baseballs (about a 
hundred of them) with a initial velocity that is too large. They are then all shot towards the wall at their varying angles. 
After they all land, velocity is decreased and the range of angles of the new set of balls that is restricted to the range of 
angles of those that made it over previously. This process is repeated iteratively forever. After a couple of minutes, a pretty
exact estimate of the optimal angle for the baseball to travel over the wall is achieved. For each set of balls launched, the
minimum of the set's angle rnage is printed on the screen. Because the range gets smaller as the rounds go on, the minimum
angle eventually approaches the optimal angle.
