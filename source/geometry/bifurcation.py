"""
:mod:`geometry` -- Generate geometry parametrically
=======================================================

    module:: geometry
    :synopsis: Generates geometries.
    moduleauthor:: Martijn Besamusca
    moduleauthor:: Ralph Erkamps
    moduleauthor:: Rick Teuthof

bifurcation.py
---------
An example geometry of a bifurcation of a vessel into two vessels.

:Authors:
    - Martijn Besamusca
    - Ralph Erkamps
    - Rick Teuthof
"""

from geometry.vessels import Vessels


def build(width, height, start_width, width1=None, width2=None, angle=40, margin=10):
    """ Builds a bifurcation.

    :param width: The width of the image and coordinate system.
    :param height: The height of the image and coordinate system.
    :param start_width: The diameter of the inlet blood vessel.
    :param width1:  The diameter of the first outlet blood vessel.
    :param width2:  The diameter of the second outlet blood vessel.
    :param angle: The angle at which the outlet vessels connect to inlet vessel.
    :param margin: The distance maintained from the wall at the top and bottom.

    :return:
        - A `Vessels` object containing the whole system.
        - A `Vessel` object containing the inlet vessel.
        - A `Vessel` object containing the first outlet vessel.
        - A `Vessel` object containing the second outlet vessel.
    """

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
