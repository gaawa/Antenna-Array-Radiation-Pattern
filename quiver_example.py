# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import coordinate_utils

poleAng = np.pi/4
elementAzimuth = np.pi/4

thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, poleAng, elementAzimuth)
phiOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth+np.pi/2)

fig = plt.figure()

x_pos = [0, 0]
y_pos = [0, 0]
z_pos = [0, 0]
x_direction = [thetaOrientVec[0], phiOrientVec[0]]
y_direction = [thetaOrientVec[1], phiOrientVec[1]]
z_direction = [thetaOrientVec[2], phiOrientVec[2]]

ax = fig.add_subplot(111, projection='3d')
ax.quiver(x_pos, y_pos, z_pos, x_direction, y_direction, z_direction)
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])
plt.show()