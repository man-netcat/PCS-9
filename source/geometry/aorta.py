# model based on https://link.springer.com/content/pdf/10.1114/1.140.pdf
# based on:
#   Normal aortoiliac diameters by CT.
#   Horejs, David; Gilbert, Perry M.; Burstein, Scott; Vogelzang, Robert L.

# Mesenteric:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3496547/

from geometry.veins import Veins


def build_thoracic(height, scale=10, margin=10):

    return


def build_abdominal(scale=10):
    """

    :param scale: px per cm
    :return:
    """

    width = scale * 15
    height = scale * 6

    w_aorta_start = 1.75
    w_aorta_end = 1.6
    pos_celiac = 1
    w_celiac = .78
    pos_lr_renal = 4
    w_lr_renal = .5
    pos_bifur = 11
    pos_superior_mesenteric = 3
    w_superior_mesenteric = .4
    pos_inferior_mesenteric = pos_superior_mesenteric + 5
    w_inferior_mesenteric = .8
    w_lr_iliac = 1.04

    pos_aorta_bot = height/2 + min(w_aorta_start, w_aorta_end)*scale/2-.3*scale
    pos_aorta_top = height/2 - min(w_aorta_start, w_aorta_end)*scale/2+.3*scale

    abdominal = Veins(width, height)
    aorta = abdominal.add_vein(
        pos_from=(0, height/2),
        pos_to=(pos_bifur*scale, height/2),
        width=w_aorta_start*scale
    )
    aorta.taper_to(w_aorta_end*scale)
    a1 = aorta.add_end(30)
    a2 = aorta.add_end(-30)

    left_iliac = aorta.append_vein((width, 5*scale), end_i=a1, width=w_lr_iliac*scale)
    right_iliac = aorta.append_vein((width, 1*scale), end_i=a2, width=w_lr_iliac*scale)

    celiac = abdominal.add_vein(
        pos_from=(pos_celiac*scale, pos_aorta_top),
        pos_to=(pos_celiac*scale + 1*scale, 0),
        width=w_celiac*scale,
        angle_from=-80,
        angle_to=-90
    )

    left_renal = abdominal.add_vein(
        pos_from=(pos_lr_renal*scale, pos_aorta_bot),
        pos_to=(pos_lr_renal*scale + 1*scale, height),
        width=w_lr_renal*scale,
        angle_from=80,
        angle_to=90)
    right_renal = abdominal.add_vein(
        pos_from=(pos_lr_renal*scale, pos_aorta_top),
        pos_to=(pos_lr_renal*scale + 1*scale, 0),
        width=w_lr_renal*scale,
        angle_from=-80,
        angle_to=-90)

    superior_mesenteric = abdominal.add_vein(
        pos_from=(pos_superior_mesenteric*scale, pos_aorta_top),
        pos_to=(pos_superior_mesenteric*scale + .5*scale, 0),
        width=w_superior_mesenteric*scale,
        angle_from=-70,
        angle_to=-90
    )
    inferior_mesenteric = abdominal.add_vein(
        pos_from=(pos_inferior_mesenteric*scale, pos_aorta_top),
        pos_to=(pos_inferior_mesenteric*scale + 1*scale, 0),
        width=w_inferior_mesenteric*scale,
        angle_from=-70,
        angle_to=-90
    )
    arteries = {
        'aorta': aorta,
        'celiac': celiac,
        'superior_mesenteric': superior_mesenteric,
        'left_renal': left_renal,
        'right_renal': right_renal,
        'inferior_mesenteric': inferior_mesenteric,
        'left_iliac': left_iliac,
        'right_iliac': right_iliac
    }
    return abdominal, arteries
