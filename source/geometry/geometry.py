"""
:mod:`geometry` -- Generate geometry parametrically
=======================================================

    module:: geometry
    :synopsis: Generates geometries.
    moduleauthor:: Martijn Besamusca
    moduleauthor:: Ralph Erkamps
    moduleauthor:: Rick Teuthof

vessel.py
---------
Low level geometry generation.
Used to generate edges and nodes.

:Authors:
    - Martijn Besamusca
    - Ralph Erkamps
    - Rick Teuthof
"""

import numpy as np


def edge(org, dest, ang_o, ang_d, width, res):
    """ Makes a smooth edge.

    :param org: Point of origin.
    :param dest: Destination point.
    :param ang_o: Angle of line at the origin.
    :param ang_d: Angle of line at the destination.
    :param width: The width of the entire edge or
                  a list of width with length res containing the widths over the entire line.
    :param res: The number of straight segments of the edge.
    :return: The coordinates of the final polygon.
    """
    org = np.array(org)
    dest = np.array(dest)
    ang_o = np.radians(ang_o)
    ang_d = np.radians(ang_d)
    l = np.linalg.norm(org - dest)
    start_w = 1 / 3 * l
    end_w = 1 / 3 * l

    p0 = org
    p1 = org + np.array([np.cos(ang_o), np.sin(ang_o)]) * start_w
    p2 = dest - np.array([np.cos(ang_d), np.sin(ang_d)]) * end_w
    p3 = dest
    if width is 0:
        return _bezier(p0, p1, p2, p3, res)
    return bezier(p0, p1, p2, p3, width, res)


def node(p, width, angle, angles, widths):
    """ Make a smooth connecting node from 1 edge to multiple edge.

    :param p: Ending point of the origin edge.
    :param width: The width of the origin edge.
    :param angle: The angle of the orign edge.
    :param angles: The angles of the connecting edges in the order of connecting (first to last).
    :param widths: The widths of the connecting edges.

    :return: The coordinates of the final polygon.
    """
    # make connecting geometry
    c_geom = [(0., 0.)]
    for a, w in zip(angles, widths):
        prev = c_geom[-1]
        a = np.radians(a - 90)
        point = (prev[0] + np.cos(a) * w, prev[1] + np.sin(a) * w)
        c_geom.append(point)
    c_geom = np.array(c_geom)

    # Move mid point of first and last to origin
    offset = c_geom[-1][0] / 2, c_geom[-1][1] / 2
    c_geom = [p - offset for p in c_geom]
    # move self width away from end
    p = np.array(p)
    perp_a = np.radians(angle + 90)
    perp = np.array([np.cos(perp_a), np.sin(perp_a)]) * width / 2
    p_mid = p + np.cross(perp, [0, 0, 1])[:2]
    c_geom = [p + p_mid for p in c_geom]

    # add ends
    p0 = p - perp
    p1 = p + perp
    c_geom.append(p0)
    c_geom.append(p1)

    xs = [p[0] for p in c_geom]
    ys = [p[1] for p in c_geom]

    return xs, ys


def bezier(p0, p1, p2, p3, width, res=10):
    """ Make a bezier curve with a given thickness.

    :param p0: Origin point.
    :param p1: Origin control point.
    :param p2: Destination control point
    :param p3: Destination point
    :param width: The width of the entire curve or
                  a list of width with length res containing the widths over the entire line.
    :param res: The number of straight segments of the edge.
    :return: The coordinates of the final polygon.
    """
    b = _bezier(p0, p1, p2, p3, res=res)  # bezier
    b_der = _bezier_deriv(p0, p1, p2, p3, res=res)  # bezier derivative
    b = np.transpose(b)
    b_der = np.transpose(b_der)

    if not hasattr(width, '__iter__'):
        width = [width] * len(b)

    c_geom = []
    for pos, der, w in zip(b, b_der, width):
        perp = np.cross(der, [0, 0, 1])[:2]
        perp = perp / np.linalg.norm(perp) * w / 2
        c_geom.append([a + b for a, b in zip(pos, perp)])
        c_geom.insert(0, [a - b for a, b in zip(pos, perp)])
    return zip(*c_geom)


def _bezier(p0, p1, p2, p3, res=10):
    """ Make a bezier curve.

    :param p0: Origin point.
    :param p1: Origin control point.
    :param p2: Destination control point
    :param p3: Destination point
    :param res: The number of straight segments of the edge.
    :return: The coordinates of the final polygon.
    """
    ts = np.linspace(0., 1., num=res)
    xs = (1 - ts) ** 3 * p0[0] + 3 * (1 - ts) ** 2 * ts * p1[0] \
         + 3 * (1 - ts) * ts ** 2 * p2[0] + ts ** 3 * p3[0]
    ys = (1 - ts) ** 3 * p0[1] + 3 * (1 - ts) ** 2 * ts * p1[1] \
         + 3 * (1 - ts) * ts ** 2 * p2[1] + ts ** 3 * p3[1]
    return xs, ys


def _bezier_deriv(p0, p1, p2, p3, res=10):
    """ The derivative of the bezier curve.

    :param p0: Origin point.
    :param p1: Origin control point.
    :param p2: Destination control point
    :param p3: Destination point
    :param res: The number of straight segments of the edge.
    :return: The coordinates of the final polygon.
    """
    ts = np.linspace(0., 1., num=res)
    xs = 3 * (1 - ts) ** 2 * (p1[0] - p0[0]) + 6 * (1 - ts) * ts * (p2[0] - p1[0]) \
         + 3 * ts ** 2 * (p3[0] - p2[0])
    ys = 3 * (1 - ts) ** 2 * (p1[1] - p0[1]) + 6 * (1 - ts) * ts * (p2[1] - p1[1]) \
         + 3 * ts ** 2 * (p3[1] - p2[1])
    return xs, ys
