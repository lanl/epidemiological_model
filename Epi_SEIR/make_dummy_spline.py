

"""NOT PART OF MAIN MODEL: FOR TESTING PURPOSES ONLY
produces a dummy spline for testing time-dependent parameter capability in
model run
"""

from numpy import linspace
from scipy.interpolate import splrep, splev
import math
import numpy as np
import pickle
import yaml

def sine(x):
    return math.sin(x) + 1.2

def main():
    with open('config/local_timedeptest_config.yaml', 'r') as in_file:
        config_dict = yaml.safe_load(in_file)

    duration = config_dict['DURATION']
    resolution = config_dict['RESOLUTION']
    print(duration, resolution)
    num_points = duration*resolution+1
    x = linspace(0, duration, num_points)
    y = np.vectorize(sine)(x) # just generates dummy data

    number_knots = int(num_points * 0.15) #number of internal knots
    knots = list(linspace(0, duration, number_knots+2))
    knots.remove(0)
    knots.remove(duration)

    spl = splrep(x, y, k=3, t=knots)
    
    # evaluate spline at time step 4.1
    #print(splev([4.1], spl)[0])

    pickle.dump(spl, open("parameters/test_spline.pkl", "wb"))
    
if __name__=="__main__":
    main()
