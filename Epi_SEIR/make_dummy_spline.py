"""NOT PART OF MAIN MODEL: FOR TESTING PURPOSES ONLY
produces a dummy spline for testing time-dependent parameter capability in
model run
"""

from numpy import linspace
from scipy.interpolate import CubicSpline
import math
import numpy as np
import pickle
import yaml

def main():
    with open('config/local_timedeptest_config.yaml', 'r') as in_file:
        config_dict = yaml.safe_load(in_file)

    duration = config_dict['DURATION']
    resolution = config_dict['RESOLUTION']
    biting_rate = 0.0014

    print(duration, resolution)

    print('x:', x)
    print('y:', y)

    spl = CubicSpline(np.linspace(0, duration, 3), biting_rate*np.ones(3))
    
    # evaluate spline at time step 4.1
    #print(spl([4.1])[0])

    pickle.dump(spl, open("parameters/test_spline.pkl", "wb"))
    
if __name__=="__main__":
    main()
