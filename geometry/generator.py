import numpy as np
import skimage.draw as draw

WALL = True
NO_WALL = False


class Veins:
    def __init__(self, res):
        self.image = []
        self.res = res
        self.veins = []

    def add_vein(self, pos_from, pos_to, width=3):
        new_vein = Vein(pos_from, pos_to, width)
        self.veins.append(new_vein)
        return new_vein

    def get_image(self):
        height = 500
        width = 1000
        image = np.full((height, width), WALL)
        for vein in self.veins:
            image = vein.draw(image)
        return image


class Vein:
    def __init__(self, pos_from, pos_to, width=5):
        self.width = width
        self.pos_from = pos_from
        self.pos_to = pos_to
        self.connects_to = []

    def append_vein(self, pos, width=5):
        new_vein = Vein(self.pos_to, pos, width)
        self.connects_to.append(new_vein)
        return new_vein

    def draw(self, image):
        for vein in self.connects_to:
            image = vein.draw(image)
        print(self.pos_from, self.pos_to)
        perp = (self.pos_from[1] - self.pos_to[1], -self.pos_from[0] + self.pos_to[0])
        print(perp)
        perp = perp / np.linalg.norm(perp) * self.width
        rect = [self.pos_from + perp, self.pos_from - perp, self.pos_to - perp, self.pos_to + perp]
        print(rect)
        rect_cc, rect_rr = zip(*rect)
        rr, cc = draw.polygon(rect_rr, rect_cc, image.shape)
        rr_circ1, cc_circ1 = draw.circle(self.pos_to[1], self.pos_to[0], self.width, image.shape)
        rr_circ2, cc_circ2 = draw.circle(self.pos_from[1], self.pos_from[0], self.width, image.shape)
        image[rr, cc] = NO_WALL
        image[rr_circ1, cc_circ1] = NO_WALL
        image[rr_circ2, cc_circ2] = NO_WALL
        return image


