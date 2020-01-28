"""
:mod:`geometry` -- Generate geometry parametrically
=======================================================

    module:: geometry
    :synopsis: Generates geometries.
    moduleauthor:: Martijn Besamusca
    moduleauthor:: Ralph Erkamps
    moduleauthor:: Rick Teuthof

vessels.py
----------
High level blood vessel geometry generation.

Basic use::
    vessels = Vessels(...)
    vessel1 = vessels.add_vessel(...)
    ...
    img = vessels.get_image()

:Authors:
    - Martijn Besamusca
    - Ralph Erkamps
    - Rick Teuthof
"""

import numpy as np
import skimage.draw as draw
from typing import List, Tuple

from geometry.geometry import bezier, edge, node

# Type indication of a position tuple
Pos = Tuple[float, float]

# Predefined values to use in the image.
WALL = True
NO_WALL = False


class Vessels:
    """ A collection of vessels.

    This collects the start of separate vessels and generates an
    image.
    """
    def __init__(self, width, height):
        """ Initialize a collection of blood vessels.

        :param width: The width of the resolution of final generated image.
        :param height: The height of the resolution of final generated image.
        """

        self.image = []
        self.width = width
        self.height = height
        self.vessels = []
        self.angles = []

    def add_vessel(self, pos_from, pos_to, angle_from=0., angle_to=0., width=3):
        """ Add a vessel to collection.

        :param pos_from: The origin point.
        :param pos_to: The destination point.
        :param angle_from: Angle at the origin.
        :param angle_to: Angle at the destination.
        :param width: The width of the entire edge or
                      a list of width with length res containing the widths over the entire line.
        :return: The new vessel.
        """

        new_vessel = Vessel(pos_from, pos_to, angle_from=angle_from, angle_to=angle_to, width=width)
        self.vessels.append(new_vessel)
        self.angles.append(angle_from)
        return new_vessel

    def get_image(self):
        """ Generates an image.

        Also updates geometry, thus this should be called before generating probe points.

        :return: An 2d array containing the image made of WALL and NO_WALL values.
        """

        image = np.full((self.height, self.width), WALL)
        for vessel in self.vessels:
            image = vessel.draw(image)
        return image


class Vessel:
    """ A single blood vessel. """

    def __init__(self, pos_from: Pos, pos_to: Pos, angle_from: float, angle_to: float, ends: List[float] = None, width: float = 5, res=50):
        """ Initialize a vessel

        Don't call this directly, use `Vessels.add_vessel` or `Vessel.append_vessel` instead.

        :param pos_from: The origin point.
        :param pos_to: The destination point.
        :param angle_from: Angle at the origin.
        :param angle_to: Angle at the destination.
        :param width: The width of the entire edge or
        :param ends: (optional) list of angles the connecting vessels have in order from start to end/
        :param width: The width of the entire edge or
                      a list of width with length res containing the widths over the entire line.
        :param res: The resolution of the edges. It is the number of segments of the bezier curve used.
        """

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

    def set_width(self, width):
        """ Set the width for the full vessel.

        :param width: New width of the entire vessel.
        """

        self.width = [width] * self.res

    def get_probe_point(self, pos):
        """ Get probe

        :param pos: Position between 0 and 1 of a probe point along the blood vessel.
        :return: A coordinate
        """

        xs, ys = edge(self.pos_from, self.pos_to, self.angle_from, self.angle_to, width=0, res=self.res)
        i = round((len(xs) - 1) * pos)
        return xs[i], ys[i]

    def taper_to(self, end_width):
        """ Tapers the blood vessel to a certain width.

        :param end_width: The final width at the end point of the vessel.
        """

        cur_width = self.start_width()
        self.width = np.linspace(cur_width, end_width, num=self.res)

    def add_narrowing(self, loc, length, scale):
        """ Narrows a part of the blood vessel.

        This is calculated via the formula: $width=width \cdot \( 2 ^ {-\(x^2\)} \)$
        This makes a good and smooth shape for a narrowing in a blood vessel.

        :param loc: A value between 0 and 1 indicating at what position the narrowing is.
        :param length: A value between 0 and 1 indicating how long the narrowing is
        :param scale: A value between 0 and 1 indicating the amount the vessel is closed.
                      0 is completely open and 1 is entirely closed.
        """

        num = round(self.res * length)
        offset = np.linspace(-2, 2, num=num)
        offset = 1 - 2 ** (-(offset ** 2)) * scale
        start = int((self.res - num) * loc)
        b = [o * w for o, w in zip(offset, self.width[start: start + num])]
        self.width[start: start + num] = b

    def add_end(self, angle):
        """ Adds a new end to the list of ends.

        This end should be used by at least one vessel, otherwise it will generate errors.
        Also this end gets added at the end of the list of ends and this angle will be drawn
        after the previously defined ends.

        :param angle: The angle of the end.
        :return: The index of the end. Is used in `append_vessel(..., end_i=...)
        """

        self._junction = None
        self.ends.append((angle, []))
        return len(self.ends) - 1

    def at_end(self, angle):
        """ Searches got an end with the given end.

        :param angle: The angle to search for.
        :return: The index of the matching end with the least number of vessels already connected.
        """

        finds = [(len(p), i) for i, (a, p) in enumerate(self.ends) if a == angle]
        finds = sorted(finds)
        return finds[0][1]

    def append_vessel(self, pos: Pos, end_i: int, angle_to=0, width: float = 5):
        """ Add a blood vessel connected to the end of this vessel.

        :param pos: The destination point.
        :param end_i: The index of the end to connect this vessel to.
        :param angle_to: Angle at the destination.
        :param width: The width of the entire edge or
                      a list of width with length res containing the widths over the entire line.
        :return: The new vessel.
        """

        self._junction = None
        new_vessel = Vessel((-1, -1), pos, angle_from=self.ends[end_i][0], angle_to=angle_to, ends=[], width=width)
        self.ends[end_i][1].append(new_vessel)
        return new_vessel

    def get_ends_loc(self, i):
        """ Get the mid point of a specific end.

        :param i: The index of the end.
        :return: The coordinate of the point.
        """

        xs, ys = self._get_junction()
        x1, x2 = xs[i:i + 2]
        y1, y2 = ys[i:i + 2]
        return (x1+x2)/2, (y1+y2)/2

    def draw(self, image):
        """ Draws this vessel and connected vessels on a image.

        :param image: An 2d array.
        :return: The image.
        """

        for i, (angle, vessels) in enumerate(self.ends):
            for vessel in vessels:
                end = self.get_ends_loc(i)
                vessel.pos_from = end
                image = vessel.draw(image)

        self._draw_vessel(image)
        self._draw_junction(image)
        return image

    def start_width(self):
        """ Get the width at the start.

        :return: The start width.
        """

        if hasattr(self.width, '__iter__'):
            return self.width[0]
        else:
            return self.width

    def end_width(self):
        """ Get the width at the end.

        :return: The end width.
        """

        if hasattr(self.width, '__iter__'):
            return self.width[-1]
        else:
            return self.width

    def _draw_vessel(self, image):
        """ Draws the edge part of the blood vessel.

        :param image: An 2d array.
        """

        xs, ys = edge(self.pos_from, self.pos_to, self.angle_from, self.angle_to, width=self.width, res=self.res)
        rr, cc = draw.polygon(ys, xs, shape=image.shape)
        image[rr, cc] = NO_WALL

    def _get_junction(self):
        """ Generates a fitting junction to connect the connected blood vessels with.

        :return: The junction polygon.
        """

        if self._junction:
            return self._junction
        angles = [end[0] for end in self.ends]
        widths = [max([v.start_width() for v in vessels]) for _, vessels in self.ends]
        xs, ys = node(self.pos_to, self.end_width(), 0, angles, widths)
        self._junction = xs, ys
        return xs, ys

    def _draw_junction(self, image):
        """ Draws the connecting node part of the blood vessel.

        :param image: An 2d array.
        """
        xs, ys = self._get_junction()
        rr, cc = draw.polygon(ys, xs, shape=image.shape)
        image[rr, cc] = NO_WALL
