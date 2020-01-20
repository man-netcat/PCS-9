import sys

def init_parameters(argv):
    
    # Fixed parameters:
    height = 80
    width = 200
    
    # Custom parameters for simulation
    viscosity = 0.05
    u0 = 0.25                           # initial and in-flow speed
    geometry = "./Geometries/geometry2.png"

    if len(sys.argv) > 1:
        viscosity = argv[1]
    if len(sys.argv) > 2:
        u0 = argv[2]
    if len(sys.argv) > 3:
        geometry = "./geometries/" + argv[3] + ".png"

    omega = 1 / (3*viscosity + 0.5)     # parameter for "relaxation"

    return (height, width, viscosity, u0, geometry, omega)