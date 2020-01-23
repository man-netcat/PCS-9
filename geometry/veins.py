import numpy as np
import skimage.draw as draw
from typing import List, Tuple

from geometry.geometry import bezier, edge, node

Pos = Tuple[float, float]

WALL = True
NO_WALL = False


class Veins:
    def __init__(self, width, height):
        self.image = []
        self.width = width
        self.height = height
        self.veins = []
        self.angles = []

    def add_vein(self, pos_from, pos_to, angle_from=0., angle_to=0., width=3):
        new_vein = Vein(pos_from, pos_to, angle_from=angle_from, angle_to=angle_to, width=width)
        self.veins.append(new_vein)
        self.angles.append(angle_from)
        return new_vein

    def get_image(self):
        image = np.full((self.height, self.width), WALL)
        for vein in self.veins:
            image = vein.draw(image)
        return image


class Vein:
    def __init__(self, pos_from: Pos, pos_to: Pos, angle_from: float, angle_to: float, ends: List[float] = None, width: float = 5, res=50):
        self.width = [width] * res
        self.pos_from = pos_from
        self.pos_to = pos_to
        self.angle_from = angle_from
        self.angle_to = angle_to
        if not ends:
            ends = []
        self.ends = [(angle, []) for angle in ends]
        self._junction = None
        self.res = res

    def get_probe_point(self, pos):
        xs, ys = edge(self.pos_from, self.pos_to, self.angle_from, self.angle_to, width=0, res=self.res)
        # pos = 1 - pos
        i = round((len(xs) - 1) * pos)
        return xs[i], ys[i]

    def add_narrowing(self, loc, width, height):
        """

        :param loc:
        :param width: (0 - 1)
        :param height: (0 - 1)
        :return:
        """
        num = round(self.res * width)
        offset = np.linspace(-2, 2, num=num)
        offset = 1 - 2 ** (-(offset ** 2)) * height
        start = int((self.res - num) * loc)
        b = [o * w for o, w in zip(offset, self.width[start: start + num])]
        self.width[start: start + num] = b

    def add_end(self, angle):
        self._junction = None
        self.ends.append((angle, []))
        return len(self.ends) - 1

    def at_end(self, angle):
        finds = [(len(p), i) for i, (a, p) in enumerate(self.ends) if a == angle]
        finds = sorted(finds)
        return finds[0][1]

    def append_vein(self, pos: Pos, end_i: int, angle_to=0, width: float = 5):
        self._junction = None
        # end = self.get_ends_loc(end_i)
        new_vein = Vein((-1, -1), pos, angle_from=self.ends[end_i][0], angle_to=angle_to, ends=[], width=width)
        self.ends[end_i][1].append(new_vein)
        return new_vein

    def get_ends_loc(self, i):
        xs, ys = self._get_junction()
        x1, x2 = xs[i:i + 2]
        y1, y2 = ys[i:i + 2]
        return (x1+x2)/2, (y1+y2)/2

    def draw(self, image):
        for i, (angle, veins) in enumerate(self.ends):
            for vein in veins:
                end = self.get_ends_loc(i)
                vein.pos_from = end
                image = vein.draw(image)

        self._draw_vein(image)
        self._draw_junction(image)
        return image

    def start_width(self):
        if hasattr(self.width, '__iter__'):
            return self.width[0]
        else:
            return self.width

    def end_width(self):
        if hasattr(self.width, '__iter__'):
            return self.width[-1]
        else:
            return self.width

    def _get_vein(self):
        return edge(self.pos_from, self.pos_to, self.angle_from, self.angle_to, width=self.width, res=self.res)

    def _draw_vein(self, image):
        xs, ys = self._get_vein()
        rr, cc = draw.polygon(ys, xs, shape=image.shape)
        image[rr, cc] = NO_WALL

    def _get_junction(self):
        if self._junction:
            return self._junction
        angles = [end[0] for end in self.ends]
        widths = [max([v.start_width() for v in veins]) for _, veins in self.ends]
        xs, ys = node(self.pos_to, self.end_width(), 0, angles, widths)
        self._junction = xs, ys
        return xs, ys

    def _draw_junction(self, image):
        xs, ys = self._get_junction()
        rr, cc = draw.polygon(ys, xs, shape=image.shape)
        image[rr, cc] = NO_WALL
