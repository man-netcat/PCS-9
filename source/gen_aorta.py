"""
gen_aorta.py
---------
Generates and saves different aorta geometries.
Use via `python3 gen_aorta.py <args>`.
For all the options type `python3 gen_aorta.py --help`.
The main option to notice is the two types:
    - Type 0 will generate geometries where just one artery is narrowed at the time
    - Type 1 will generate all combinations of narrowings.

:Authors:
    - Martijn Besamusca
    - Ralph Erkamps
    - Rick Teuthof
"""


import argparse
import os
import numpy as np
import shutil
from PIL import Image
import itertools
import matplotlib.pyplot as plt

import geometry.aorta as aorta

def parse_args():
    """ Parse given arguments.

    :return: Dictionary containing all passed options.
    """
    # parse arguments
    parser = argparse.ArgumentParser(description='Generate multiple bifurcation geom files')
    parser.add_argument('--type', '-t', type=int, default=0, help='0 for just one artery change width at the time'
                                                                  '1 for all combinations of widths')
    parser.add_argument('--res', '-R', type=int, default=20, help='Number of pixels per centimeter')
    parser.add_argument('--num', '-n', type=int, default=3, help='Number of steps per vein')
    parser.add_argument('--min', '-m', type=float, default=0.0, help='Min width of the narrowing')
    parser.add_argument('--max', '-M', type=float, default=0.5, help='Max width of the narrowing')
    parser.add_argument('--out', '-o', type=str, default='out', help='The name of the output folder')
    parser.add_argument('--clear', '-c', action='store_true', help='Clear the output folder')
    parser.add_argument('--draw', '-d', action='store_true', help='Draw probes at end')
    return parser.parse_args()


def make_output_folder(args):
    """ Make the output directory if it doesn't exist.

    :param args: Passed options.
    """
    # make folder
    if not args.out:
        print('No out path')
        exit(1)
    folder = args.out
    if args.clear and os.path.exists(folder):
        shutil.rmtree(folder)
        os.mkdir(folder)
    elif not os.path.exists(folder):
        os.mkdir(folder)


def generate(args):
    """ Generate all geometries.

    :param args: Passed options.
    :return:
        - The Vessels object containing collection of arteries
        - A dictionary containing Vessel objects.
    """
    abdominal, arteries = aorta.build_abdominal(args.res)
    iterator = comb(args, arteries) if args.type is 1 else prod(args, arteries)
    for filename in iterator:
        print(filename)
        image = abdominal.get_image()
        new_im = Image.fromarray(image)
        new_im.save(os.path.join(args.out, filename + '.png'))
    return abdominal, arteries


def comb(args, arteries):
    """ Generate type 1 options geometries

    :param args: Passed options.
    :param arteries: A dictionary containing Vessel objects.
    :return: Yields the name of this geometry configuration.
    """
    combinations = itertools.product(range(args.num), repeat=len(arteries))
    artery_names = list(arteries.keys())
    artery_names.remove('aorta')
    artery_widths_org = {key: val.start_width() for key, val in arteries.items()}
    for widths in combinations:
        artery_widths = {name: args.min + (args.max - args.min)*(w/(args.num - 1)) for name, w in zip(artery_names, widths)}
        # yield artery_widths
        for name, width in artery_widths.items():
            if name is 'aorta' and True:
                continue
            artery = arteries.get(name)
            artery.set_width(artery_widths_org.get(name))
            artery.add_narrowing(loc=0.5, length=0.3, scale=width)
        yield '_'.join([f'{w:.4f}' for w in artery_widths.values()])


def prod(args, arteries):
    """ Generate type 0 options geometries

    :param args: Passed options.
    :param arteries: A dictionary containing Vessel objects.
    :return: Yields the name of this geometry configuration.
    """
    artery_widths_org = {key: val.start_width() for key, val in arteries.items()}
    for artery_name, artery in arteries.items():
        if artery_name is 'aorta' and True:
            continue
        for width in np.linspace(args.min, args.max, num=args.num):
            artery.add_narrowing(loc=0.5, length=0.3, scale=width)
            # yield artery_name, width
            yield f'{artery_name}_{width:.4f}'
            artery.set_width(artery_widths_org.get(artery_name))


def make_probes(args, arteries):
    """ Places all measure points.

    :param args: Passed options.
    :param arteries: A dictionary containing Vessel objects.
    :return: None
    """
    names = {
        'supraceliac_aorta': 'supraceliac aorta',
        'aorta': 'abdominal aorta',
        'celiac': 'celiac artery',
        'superior_mesenteric': 'superior mesenteric artery',
        'left_renal': 'left renal artery',
        'right_renal': 'right renal artery',
        'inferior_mesenteric': 'inferior mesenteric artery',
        'left_iliac': 'left iliac artery',
        'right_iliac': 'right iliac artery'
    }
    start_suffix = ' start'
    end_suffix = ' end'

    ax=None
    if args.draw:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        image = abdominal.get_image()
        ax.imshow(image, cmap=plt.cm.gray)

    with open(os.path.join(args.out, 'probe.txt'), 'w') as f:

        for name, artery in arteries.items():
            name = names.get(name)
            if name == 'supraceliac aorta':
                write_probe(f, artery, [0.1], [name], ax)
            elif name == 'abdominal aorta':
                write_probe(f, artery, [1], [name], ax)
            else:
                write_probe(f, artery, [0.2, 0.8], [name+start_suffix, name+end_suffix], ax)

    if args.draw:
        plt.show()


def write_probe(file, artery, locs, names, ax=None):
    """ Write a probe location to a file (and to a plot)

    :param file: A open writable file
    :param artery: The artery which is getting measured.
    :param locs: Locations measured.
    :param names: Names of the locations measured.
    :param ax: (optional) The axes to draw to.
    :return: None
    """
    probes = [[name, *artery.get_probe_point(loc)] for name, loc in zip(names, locs)]
    for name, x, y in probes:
        file.write(f'{name},{x},{y}\n')

    if ax:
        print(probes)
        print([*zip(*probes)][1:])
        labels, xs, ys = [*zip(*probes)]
        ax.scatter(xs, ys, s=20)
        for name, x, y in probes:
            ax.annotate(name, (x, y), c='C0')


if __name__ == '__main__':
    args = parse_args()
    make_output_folder(args)
    abdominal, arteries = generate(args)
    make_probes(args, arteries)
