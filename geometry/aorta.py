# model based on https://link.springer.com/content/pdf/10.1114/1.140.pdf
# based on:
#   Normal aortoiliac diameters by CT.
#   Horejs, David; Gilbert, Perry M.; Burstein, Scott; Vogelzang, Robert L.

from geometry.veins import Veins


def build_thoracic(height, scale=10, margin=10):

    return


def build_abdominal(scale=10):
    """

    :param scale: px per cm
    :return:
    """

    width = scale * 20
    height = scale * 6
    abdominal = Veins(width, height)
    aorta = abdominal.add_vein((0, height/2), (15*scale, height/2), width=1.75*scale)
    # aorta.taper_to(1.60*scale)
    a1 = aorta.add_end(30)
    a2 = aorta.add_end(-20)
    left_iliac = aorta.append_vein((width, 5*scale), end_i=a1, width=1.04*scale)
    right_iliac = aorta.append_vein((width, 1*scale), end_i=a2, width=1.04*scale)
    celiac = abdominal.add_vein((2*scale, height/2), (2*scale, height), width=0.8*scale, angle_from=-90, angle_to=90)
    # superior_mesenteric = abdominal.add_vein((0, height/2), (10*scale, height/2), width=2.54*scale)
    # left_renal = abdominal.add_vein((0, height/2), (10*scale, height/2), width=2.54*scale)
    # right_renal = abdominal.add_vein((0, height/2), (10*scale, height/2), width=2.54*scale)
    # inferior_mesenteric = abdominal.add_vein((0, height/2), (10*scale, height/2), width=2.54*scale)
    return abdominal