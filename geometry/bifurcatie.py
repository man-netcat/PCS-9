from geometry.veins import Veins


def build(width, height, start_width, width1=None, width2=None, angle=40):
    # if width1 is None: width1 = width / 2
    bifurcation = Veins(width, height)
    start = (0, height / 2)
    mid = (width / 3, height / 2)
    ends = [(width, height / 3), (width, height / 3 * 2)]
    v0 = bifurcation.add_vein(start, mid, width=start_width)
    a1 = v0.add_end(angle)
    a0 = v0.add_end(-angle)
    v0.append_vein(ends[0], end_i=a0, width=start_width )
    v0.append_vein(ends[1], end_i=a1, width=start_width )
    return bifurcation
