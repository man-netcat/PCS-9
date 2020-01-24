import argparse
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
from PIL import Image
import itertools

import geometry.aorta as aorta

def parse_args():
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
    abdominal, arteries = aorta.build_abdominal(args.res)
    iterator = comb(args, arteries) if args.type is 1 else prod(args, arteries)
    for filename in iterator:
        print(filename)
        image = abdominal.get_image()
        new_im = Image.fromarray(image)
        new_im.save(os.path.join(args.out, filename + '.png'))


def comb(args, arteries):
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
            artery.add_narrowing(loc=0.5, width=0.3, height=width)
        yield '_'.join([f'{w:.4f}' for w in artery_widths.values()])


def prod(args, arteries):
    artery_widths_org = {key: val.start_width() for key, val in arteries.items()}
    for artery_name, artery in arteries.items():
        if artery_name is 'aorta' and True:
            continue
        for width in np.linspace(args.min, args.max, num=args.num):
            artery.add_narrowing(loc=0.5, width=0.3, height=width)
            # yield artery_name, width
            yield f'{artery_name}_{width:.4f}'
            artery.set_width(artery_widths_org.get(artery_name))


if __name__ == '__main__':
    args = parse_args()
    make_output_folder(args)
    generate(args)



# for width in widths:
#     split_vein, arteries = aorta.build_abdominal(args.res)
#     v1.add_narrowing(loc=0.5, width=0.4, height=width)
#     image = split_vein.get_image()
#     filename = f'bifurcation_{width:.4f}.png'
#     new_im = Image.fromarray(image)
#     new_im.save(os.path.join(folder, filename))
#
# if not arteries:
#     print('Nothing generated')
#     exit(1)

# # display
# image = split_vein.get_image()
#
# probe_in = v0.get_probe_point(0.3)
# probe_normal = v2.get_probe_point(0.8)
# probe_nar_before = v1.get_probe_point(0.2)
# probe_nar_after = v1.get_probe_point(0.8)
#
# with open(os.path.join(folder, 'probe.txt'), 'w') as f:
#     s = '{},{},{}\n'
#     f.write(s.format('inlet', *probe_in))
#     f.write(s.format('normal vein', *probe_normal))
#     f.write(s.format('start narrow vein', *probe_nar_before))
#     f.write(s.format('end narrow vein', *probe_nar_after))
#     # for x, y in probe_in, probe_normal, probe_nar_before, probe_nar_after:
#     #     f.write(f'{x} {y}\n')
#
# # print(probe_in, probe_normal, probe_nar_before, probe_nar_after)
# if args.draw:
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.imshow(image, cmap=plt.cm.gray)
#     ax.scatter(*zip(*[probe_in, probe_normal, probe_nar_before, probe_nar_after]), s=20)
#     plt.show()
