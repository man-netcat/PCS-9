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

    def add_vein(self, pos_from, pos_to, angle=0., width=3):
        new_vein = Vein(pos_from, pos_to, angle=angle, width=width)
        self.veins.append(new_vein)
        self.angles.append(angle)
        return new_vein

    def get_image(self):
        image = np.full((self.height, self.width), WALL)
        for vein in self.veins:
            image = vein.draw(image)
        return image


class Vein:
    def __init__(self, pos_from: Pos, pos_to: Pos, angle: float, ends: List[float] = None, width: float = 5):
        self.width = width
        self.pos_from = pos_from
        self.pos_to = pos_to
        self.angle = angle
        if not ends:
            ends = []
        self.ends = [(angle, []) for angle in ends]
        self._junction = None

    def add_end(self, angle):
        self._junction = None
        self.ends.append((angle, []))
        return len(self.ends) - 1

    def at_end(self, angle):
        finds = [(len(p), i) for i, (a, p) in enumerate(self.ends) if a == angle]
        finds = sorted(finds)
        return finds[0][1]

    def append_vein(self, pos: Pos, end_i: int, width: float = 5):
        self._junction = None
        # end = self.get_ends_loc(end_i)
        new_vein = Vein((-1, -1), pos, self.ends[end_i][0], [], width)
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

    def _draw_vein(self, image):
        xs, ys = edge(self.pos_from, self.pos_to, self.angle, 0, width=self.width)
        rr, cc = draw.polygon(ys, xs, shape=image.shape)
        image[rr, cc] = NO_WALL

    def _get_junction(self):
        if self._junction:
            return self._junction
        angles = [end[0] for end in self.ends]
        widths = [max([v.width for v in veins]) for _, veins in self.ends]
        xs, ys = node(self.pos_to, self.width, 0, angles, widths)
        self._junction = xs, ys
        return xs, ys

    def _draw_junction(self, image):
        xs, ys = self._get_junction()
        rr, cc = draw.polygon(ys, xs, shape=image.shape)
        image[rr, cc] = NO_WALL
