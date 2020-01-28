from geometry.vessels import Vessels


def build(width, height, start_width, width1=None, width2=None, angle=40, margin=10):
    if width1 is None: width1 = start_width / 2
    if width2 is None: width2 = start_width / 2
    bifurcation = Vessels(width, height)
    start = (0, height / 2)
    mid = (width / 3, height / 2)
    ends = [(width, margin+width1), (width, height - margin - width2)]
    v0 = bifurcation.add_vessel(start, mid, width=start_width)
    a1 = v0.add_end(angle)
    a0 = v0.add_end(-angle)
    v1 = v0.append_vessel(ends[0], end_i=a0, width=width1)
    v2 = v0.append_vessel(ends[1], end_i=a1, width=width2)
    return bifurcation, v0, v1, v2
