import sys

def init_parameters(argv):
    
    # Fixed parameters:
    height = 80
    width = 200
    
    # Custom parameters for simulation
    viscosity = 0.05
    u0 = 0.25                           # initial and in-flow speed
    geometry = "./geometries/geometry2.png"
    filename = "results.txt"
    if len(sys.argv) > 1:
        geometry = "./geometries/" + argv[1] + ".png"
    if len(sys.argv) > 2:
        filename = argv[2]
    if len(sys.argv) > 3:
        u0 = argv[3]
    if len(sys.argv) > 4:
        viscosity = argv[4]
        

    omega = 1 / (3*viscosity + 0.5)     # parameter for "relaxation"

    return (height, width, viscosity, u0, geometry, omega)