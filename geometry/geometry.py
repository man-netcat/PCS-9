import numpy as np


def edge(org, dest, ang_o, ang_d, width):
    org = np.array(org)
    dest = np.array(dest)
    ang_o = np.radians(ang_o)
    ang_d = np.radians(ang_d)
    l = np.linalg.norm(org - dest)
    start_w = 1/3 * l
    end_w = 1/3 * l

    p1 = org
    p2 = org + np.array([np.cos(ang_o), np.sin(ang_o)]) * start_w
    p3 = dest - np.array([np.cos(ang_d), np.sin(ang_d)]) * end_w
    p4 = dest

    return bezier(p1, p2, p3, p4, width)
    # perp = (p1[1] - p2[1], p2[0] - p1[0])
    # perp = perp / np.linalg.norm(perp) * width
    # rect = [p1 + perp, p1 - perp, p2 - perp, p2 + perp]
    # return rect


def node(p, width, angle, angles, widths):
    # make connecting geometry
    c_geom = [(0, 0)]
    for angle, width in zip(angles, widths):
        prev = c_geom[-1]
        angle = np.radians(angle)
        point = (prev[0] + np.cos(angle) * width, prev[1] + np.sin(angle) * width)
        c_geom.append(point)
    c_geom = np.array(c_geom)
    print(c_geom)

    # Move mid point of first and last to origin
    offset = c_geom[-1][0] / 2, c_geom[-1][1] / 2
    c_geom = [p - offset for p in c_geom]
    # move self width away from end
    p = np.array(p)
    perp_a = np.radians(angle+90)
    perp = np.array([np.cos(perp_a), np.sin(perp_a)]) * width
    p_mid = p + np.cross(perp, [0, 0, 1])[:2]
    c_geom = [p + p_mid for p in c_geom]

    # add ends
    p1 = p - perp
    p2 = p + perp
    c_geom = []
    c_geom.append(p1)
    c_geom.append(p2)
    c_geom.append(p_mid)
    print(p1, p2)

    xs = [p[0] for p in c_geom]
    ys = [p[1] for p in c_geom]

    print(c_geom)
    return xs, ys
pass
    # perp = (p1[1] - p2[1], p2[0] - p1[0])
    # perp = perp / np.linalg.norm(perp) * width
    # rect = [p1 + perp, p1 - perp, p2 - perp, p2 + perp]
    # return rect


def bezier(p1, p2, p3, p4, width, res=10):
    perp = (p1[1] - p2[1], p2[0] - p1[0])
    perp = perp / np.linalg.norm(perp) * width / 2

    p1u, p2u, p3u, p4u = [p + perp for p in [p1, p2, p3, p4]]
    xsu, ysu = _bezier(p1u, p2u, p3u, p4u, res=res)
    xsd = [p1[0] - perp[0]]
    ysd = [p1[1] - perp[1]]
    for i, (x, y) in enumerate(zip(xsu[1:-1], ysu[1:-1])):
        x_old, y_old = xsu[i], ysu[i]
        p = (y - y_old, x_old - x)
        p = p / np.linalg.norm(p) * width
        xsd.append(x + p[0])
        ysd.append(y + p[1])
    xsd.append(p4[0] - perp[0])
    ysd.append(p4[1] - perp[1])

    xsd.reverse()
    ysd.reverse()

    xs = np.concatenate((xsu, xsd))
    ys = np.concatenate((ysu, ysd))

    print(xs)
    print(ys)
    return xs, ys


def _bezier(p1, p2, p3, p4, res=10):
    ts = np.linspace(0., 1., num=res)
    xs = (1-ts)**3 * p1[0] + 3*(1-ts)**2 * ts * p2[0] \
         + 3*(1-ts) * ts**2 * p3[0] + ts**3 * p4[0]
    ys = (1-ts)**3 * p1[1] + 3*(1-ts)**2 * ts * p2[1] \
         + 3*(1-ts) * ts**2 * p3[1] + ts**3 * p4[1]
    return xs, ys

def draw(image, polygons, val):
    pass
