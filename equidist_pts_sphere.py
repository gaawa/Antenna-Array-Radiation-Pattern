# source:
# https://www.youtube.com/watch?v=rwV_dAlWWnw
# https://scholar.rose-hulman.edu/rhumj/vol18/iss2/5/

import numpy as np

def spherical_coordinate(azimuth, elevation):
    xCoord = np.cos(azimuth)*np.cos(elevation)
    yCoord = np.sin(azimuth)*np.cos(elevation)
    zCoord = np.sin(elevation)
    return np.array([xCoord, yCoord, zCoord])

def generate_points(Npts):
    x_val = 0.1 + 1.2 * Npts # see the soruce paper
    pts = []
    s_start = -1 + 1/(Npts-1)
    s_steps = (2 - 2/(Npts-1))/(Npts-1)
    
    for i in range(Npts):
        s = s_start + i*s_steps
        azimuth = s*x_val
        elevation = np.pi/2 * np.sign(s) * (1-np.sqrt(1-np.abs(s)))
        pts.append(spherical_coordinate(azimuth, elevation))
        
    return pts

if __name__ == "__main__":
    # show 3d scatter plot of points as a demo
    import matplotlib.pyplot as plt
    pts = generate_points(1000)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for p in pts:
        ax.scatter(p[0], p[1], p[2], marker='o', color='b')
        
    plt.show()
        