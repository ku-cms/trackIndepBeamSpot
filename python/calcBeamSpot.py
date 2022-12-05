# calcBeamSpot.py
import numpy as np

# given parameters (b, r_2D, r_3D), calculate beam spot (x, y, z)
def calcBeamSpot(params):
    b    = params[0]
    r_2D = params[1]
    r_3D = params[2]
    phi = b + np.pi/2
    x = r_2D * np.cos(phi)
    y = r_2D * np.sin(phi)
    z = np.sqrt(r_3D**2 - r_2D**2) 
    return [x, y, z]

# given beam spot (x, y, z), calculate parameters (b, r_2D, r_3D)
def calcParams(beam_spot):
    x = beam_spot[0]
    y = beam_spot[1]
    z = beam_spot[2]
    # avoid dividing by 0
    if x != 0:
        phi = np.arctan(y / x)
        b   = phi - np.pi/2
        # stay within [-pi, pi]
        if b < -np.pi:
            b += 2*np.pi
    else:
        print("ERROR: Cannot divide by x = {0}".format(x))
        b = -999
    r_2D = np.sqrt(x**2 + y**2)
    r_3D = np.sqrt(x**2 + y**2 + z**2)
    return [b, r_2D, r_3D]

# print beam spot (x, y, z)
def printBeamSpot(beam_spot):
    x = beam_spot[0]
    y = beam_spot[1]
    z = beam_spot[2]
    print("x = {0:.3f}".format(x))
    print("y = {0:.3f}".format(y))
    print("z = {0:.3f}".format(z))

# print parameters (b, r_2D, r_3D)
def printParams(params):
    b    = params[0]
    r_2D = params[1]
    r_3D = params[2]
    print("b    = {0:.3f}".format(b))
    print("r_2D = {0:.3f}".format(r_2D))
    print("r_3D = {0:.3f}".format(r_3D))

def run(dataset):
    name      = dataset["name"]
    beam_spot = dataset["beam_spot"]
    params    = calcParams(beam_spot)
    print(name)
    printBeamSpot(beam_spot)
    printParams(params)

def main():
    # 2021 express data, run 346512
    #x =  0.1727
    #y = -0.1913
    #z =  0.3813
    #beam_spot       = [x, y, z]
    #params          = calcParams(beam_spot)
    #beam_spot_v2    = calcBeamSpot(params)
    #printBeamSpot(beam_spot)
    #printParams(params)
    #printBeamSpot(beam_spot_v2)

    datasets = []

    # TTBar MC (z-smeared)
    #x =  0.10
    #y = -0.08
    #z =  0.00
    #beam_spot = [x, y, z]
    #dataset = {}
    #dataset["name"] = "TTBar MC (z-smeared)"
    #dataset["beam_spot"] = beam_spot
    #datasets.append(dataset)

    # Single Muon 2017B
    #x =  0.09
    #y = -0.04
    #z =  0.15
    #beam_spot = [x, y, z]
    #dataset = {}
    #dataset["name"] = "Single Muon 2017B"
    #dataset["beam_spot"] = beam_spot
    #datasets.append(dataset)
    
    # Single Muon 2017B (PvZ < 0)
    x =  0.085
    y = -0.036
    z =  0.147
    beam_spot = [x, y, z]
    dataset = {}
    dataset["name"] = "Single Muon 2017B (PvZ < 0)"
    dataset["beam_spot"] = beam_spot
    datasets.append(dataset)
    
    # Single Muon 2017B (PvZ >= 0)
    x =  0.086
    y = -0.036
    z =  0.147
    beam_spot = [x, y, z]
    dataset = {}
    dataset["name"] = "Single Muon 2017B (PvZ >= 0)"
    dataset["beam_spot"] = beam_spot
    datasets.append(dataset)

    # Zero Bias 2017B
    #x =  0.09
    #y = -0.04
    #z =  0.17
    #beam_spot = [x, y, z]
    #dataset = {}
    #dataset["name"] = "Zero Bias 2017B"
    #dataset["beam_spot"] = beam_spot
    #datasets.append(dataset)

    # Express Data 2021, run 346512
    #x =  0.1727
    #y = -0.1913
    #z =  0.3813
    #x =  0.17
    #y = -0.19
    #z =  0.38
    #beam_spot = [x, y, z]
    #dataset = {}
    #dataset["name"] = "Express Data 2021"
    #dataset["beam_spot"] = beam_spot
    #datasets.append(dataset)

    for dataset in datasets:
        run(dataset)

if __name__ == "__main__":
    main()

