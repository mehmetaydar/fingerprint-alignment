"""
apply
  gabor : images/ppf-test.png 16 --save
  thinning: ppf-test_enhanced.gif --save
  extract minutias: ppf-test_enhanced_thinned.gif --save
  calculate poincare: ppf-test_enhanced_thinned_minutiae.gif 16 1 --smooth --save
  align minutias with poincare index
"""
from PIL import Image, ImageDraw
import utils
import argparse
import math
from math import sin, cos, pi, atan2
import frequency
import os
from gabor import gabor_kernel
from thining import reverse, apply_all_structures
import pickle
import json
import os.path

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt0

signum = lambda x: -1 if x < 0 else 1
cells = [(-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]

def gabor(im, W, angles):
    (x, y) = im.size
    im_load = im.load()

    freqs = frequency.freq(im, W, angles)
    print "computing local ridge frequency done"

    gauss = utils.gauss_kernel(3)
    utils.apply_kernel(freqs, gauss)

    for i in range(1, x / W - 1):
        for j in range(1, y / W - 1):
            kernel = gabor_kernel(W, angles[i][j], freqs[i][j])
            for k in range(0, W):
                for l in range(0, W):
                    im_load[i * W + k, j * W + l] = utils.apply_kernel_at(
                        lambda x, y: im_load[x, y],
                        kernel,
                        i * W + k,
                        j * W + l)

    return im

def make_thin(im):
    loaded = utils.load_image(im)
    utils.apply_to_each_pixel(loaded, lambda x: 0.0 if x > 10 else 1.0)
    print "loading phase done"

    t1 = [[1, 1, 1], [0, 1, 0], [0.1, 0.1, 0.1]]
    t2 = utils.transpose(t1)
    t3 = reverse(t1)
    t4 = utils.transpose(t3)
    t5 = [[0, 1, 0], [0.1, 1, 1], [0.1, 0.1, 0]]
    t7 = utils.transpose(t5)
    t6 = reverse(t7)
    t8 = reverse(t5)

    thinners = [t1, t2, t3, t4, t5, t6, t7]

    usage = True
    while(usage):
        usage = apply_all_structures(loaded, thinners)
        print "single thining phase done"

    print "thining done"

    utils.apply_to_each_pixel(loaded, lambda x: 255.0 * (1 - x))
    utils.load_pixels(im, loaded)
    return im

def minutiae_at(pixels, i, j):
    values = [pixels[i + k][j + l] for k, l in cells]

    crossings = 0
    for k in range(0, 8):
        crossings += abs(values[k] - values[k + 1])
    crossings /= 2

    if pixels[i][j] == 1:
        if crossings == 1:
            return "ending"
        if crossings == 3:
            return "bifurcation"
    return "none"
def calculate_minutiaes(im):
    minutiaes = []
    pixels = utils.load_image(im)
    utils.apply_to_each_pixel(pixels, lambda x: 0.0 if x > 10 else 1.0)

    (x, y) = im.size
    result = im.convert("RGB")

    draw = ImageDraw.Draw(result)

    colors = {"ending" : (150, 0, 0), "bifurcation" : (0, 150, 0)}

    indexes = {"ending" : 1, "bifurcation" : 2}

    ellipse_size = 2
    for i in range(1, x - 1):
        for j in range(1, y - 1):
            minutiae = minutiae_at(pixels, i, j)
            if minutiae != "none":
                minutiaes.append( (i, j, indexes[minutiae]) )
                draw.ellipse([(i - ellipse_size, j - ellipse_size), (i + ellipse_size, j + ellipse_size)], outline = colors[minutiae])

    del draw

    return result, minutiaes

def get_angle(left, right):
    angle = left - right
    if abs(angle) > 180:
        angle = -1 * signum(angle) * (360 - abs(angle))
    return angle

def poincare_index_at(i, j, angles, tolerance):
    deg_angles = [math.degrees(angles[i - k][j - l]) % 180 for k, l in cells]
    index = 0
    for k in range(0, 8):
        if abs(get_angle(deg_angles[k], deg_angles[k + 1])) > 90:
            deg_angles[k + 1] += 180
        index += get_angle(deg_angles[k], deg_angles[k + 1])

    if 180 - tolerance <= index and index <= 180 + tolerance:
        return "loop"
    if -180 - tolerance <= index and index <= -180 + tolerance:
        return "delta"
    if 360 - tolerance <= index and index <= 360 + tolerance:
        return "whorl"
    return "none"

def calculate_singularities(im, angles, tolerance, W, singularity_type = None):
    cores = []
    (x, y) = im.size
    result = im.convert("RGB")

    draw = ImageDraw.Draw(result)

    #loop = red, delta = green, whorl = blue
    colors = {"loop" : (150, 0, 0), "delta" : (0, 150, 0), "whorl": (0, 0, 150)}

    center = (len(angles)/2, len(angles)/2)
    for i in range(1, len(angles) - 1):
        for j in range(1, len(angles[i]) - 1):
            singularity = poincare_index_at(i, j, angles, tolerance)
            if singularity != "none":
                if singularity_type is not None and singularity != singularity_type:
                    continue
                distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip((i, j), center)]))
                cores.append((distance, i * W, j * W, singularity) )
                draw.ellipse([(i * W, j * W), ((i + 1) * W, (j + 1) * W)], outline = colors[singularity])

    cores = sorted(cores)
    del draw

    return result, cores

def chooseCore(cores, bsize):
    if len(cores) == 0: return None

    core = cores[0]
    return atan2(core[2], core[1]), core[1], core[2], core[3]

    sels = {}
    for index in range(1, len(cores)-1):
        core0 = cores[index-1]
        core1 = cores[index]

        i0 = core0[1]/bsize
        j0 = core0[2]/bsize
        i1 = core1[1]/bsize
        j1 = core1[2]/bsize

        if core1[3] == core0[3] and abs(i1-i0)<=1 and abs(j1-j0)<=1:
            sels[i0,j0] = core0
            sels[i1, j1] = core1
        else:
            break
    if len(sels) > 0:
        ib=0
        jb=0
        for key, core in sels.iteritems():
            ib+= core[1]
            jb += core[2]
        y = ib/len(sels)
        x = jb/len(sels)
        return atan2(y,x), y, x, core[3]
    else:
        core = cores[0]
        return atan2(core[2], core[1]), core[1], core[2], core[3]

def alignMinutias(core, minutiaes):
    if core is None:
        print "No core to align!"
        return minutiaes
    dir = core[0] #degree in radiant
    corex = core[1]
    corey = core[2]
    radcorexy = atan2(corey, corex)
    minutiaes2 = []
    for minutia in minutiaes:
        #x = minutia[0]
        #y = minutia[1]
        x = minutia[0]-corex
        y = minutia[1]-corey
        xy = np.array([x, y])
        radxy = atan2(y, x)
        #x2 = x - corex
        #y2 = y - corey
        #x3 = x2 * cos(dir) + y2 * sin(dir)
        #y3 = -x2 * sin(dir) + y2 * cos(dir)

        #th = -radxy #counter clockwise
        #th = radxy
        th = radcorexy
        c, s = cos(th), sin(th)
        R = np.array(((c, -s), (s, c)))
        #newxy = np.matmul(R, xy).tolist()
        newxy = R.dot(xy).tolist()
        #print newxy
        minutiaes2.append( ( int(round(newxy[0])), int(round(newxy[1])), minutia[2]) )
    return minutiaes2

#not being used
def alignBounds(core, bounds):
    if core is None:
        print "No core to align!"
        return bounds
    (minx, maxx, miny, maxy) = bounds
    points = [(0, 0), (0, maxy), (maxx, maxy), (maxx, 0)]
    dir = core[0] #degree in radiant
    corex = core[1]
    corey = core[2]
    points2 = []
    for point in points:
        x = point[0]
        y = point[1]
        radxy = atan2(y, x)
        x2 = x - corex
        y2 = y - corey
        x3 = x2 * cos(dir) + y2 * sin(dir)
        y3 = -x2 * sin(dir) + y2 * cos(dir)
        points2.append( ( int(round(x3)), int(round(y3))) )

    minx2 = min([xx[0] for xx in points2])
    maxx2 = max([xx[0] for xx in points2])
    miny2 = min([xx[1] for xx in points2])
    maxy2 = max([xx[1] for xx in points2])
    bounds2 = (minx2, maxx2, miny2, maxy2)
    return bounds2

#not being used
def matrixMultiply(X, Y):
    result = [[sum(a * b for a, b in zip(X_row, Y_col)) for Y_col in zip(*Y)] for X_row in X]
    return result
def plotAndSaveMinutiae(minutiaes, name):
    print "Plotting and saving minutiaes: ", name
    mins0 = [[m[0], m[1]] for m in minutiaes]
    data = np.array([
        mins0
    ])
    x, y = data.T
    plt0.scatter(x, y, s=1, c='black')
    #plt0.axis('off')
    plt0.savefig(name)
    plt0.close()

def plotAndSaveMinutiasPIL(minutiaes, name, size_xy):
    sizex = size_xy[0]
    sizey = size_xy[1]
    xx = [m[0] for m in minutiaes]
    yy = [m[1] for m in minutiaes]
    minx = abs(min(xx))+5
    miny = abs(min(yy))+5
    newminutiaes = [(m[0]+abs(minx), m[1]+abs(miny), m[2]) for m in minutiaes]
    newsize_xy = ( sizex+abs(minx), sizey+abs(miny) )
    img = Image.new('RGB', newsize_xy, color='white')
    #result = img.convert("RGB")

    draw = ImageDraw.Draw(img)

    colors = {1: (150, 0, 0), 2: (0, 150, 0)}
    indexes = {"ending": 1, "bifurcation": 2}

    ellipse_size = 1
    for minutiae in newminutiaes:
        i = minutiae[0]
        j = minutiae[1]
        col = minutiae[2]
        outlinec = colors[col]
        draw.ellipse([(i - ellipse_size, j - ellipse_size), (i + ellipse_size, j + ellipse_size)],
                     outline=outlinec)
    img.show()
    img.save(name)
    del draw

def main1():
    # alignMinutias((2,8, 'loop'), [[40, 120, 'bip']])
    parser = argparse.ArgumentParser(
        description="Gabor, thinning, extract minutias, calculate poincare and align minutias.")
    parser.add_argument("image", nargs=1, help="Path to image")
    # parser.add_argument("block_size", nargs=1, help = "Block size")
    parser.add_argument("tolerance", nargs=1, help="Tolerance for Poincare index")
    parser.add_argument("block_size", nargs=1, help="Blocksize")
    parser.add_argument('--preprocess', "-p", action='store_true',
                        help="Preprocess the image with: Gabor filtering and thinning")
    parser.add_argument('--smooth', "-s", action='store_true', help="Use Gauss for smoothing")
    parser.add_argument("--save", action='store_true', help="Save result image as src_image_enhanced.gif")
    args = parser.parse_args()

    #im = Image.open(args.image[0])
    # im = Image.open("ppf-test_enhanced_thinned.gif")
    #im = im.convert("L")  # covert to grayscale
    # im.show()

    W = int(args.block_size[0])
    #W = 16
    print "Block-size: ", W
    singularity_type = None

    f = lambda x, y: 2 * x * y
    g = lambda x, y: x ** 2 - y ** 2

    files = []
    for i in range(1, 4):
        for j in range(1, 4):
            if i < 10:
                finputminutia = '/Users/Admin/fvs/samples/DB3_B/10' + str(i) + '_' + str(j) + '.tif'
            else:
                finputminutia = '/Users/Admin/fvs/samples/DB3_B/1' + str(i) + '_' + str(j) + '.tif'
            files.append({
                'i': i,
                'j': j,
                'image': finputminutia
            })
    ffirst_singularity = {}
    print "Gabor filter, thinning, pointcare, extract minutias, alignment, exporting minutias to file ..."
    for file in files:
        i = file['i']
        j = file['j']
        print "image: ", os.path.splitext(file['image'])
        base_image_name = os.path.splitext(os.path.basename(file['image']))[0]
        image_enhanced_loc = "/Users/Admin/fvs/samples/dumps2/" + base_image_name + "_enhanced.gif"
        image_enhanced_minutiae_loc = "/Users/Admin/fvs/samples/dumps2/" + base_image_name + "_enhanced_minutiaes.png"
        image_enhanced_minutiae_aligned_loc = "/Users/Admin/fvs/samples/dumps2/" + base_image_name + "_enhanced_minutiaes_aligned.png"
        print "image enhanced: ", image_enhanced_loc
        # f['image']
        if args.preprocess or os.path.isfile(image_enhanced_loc) == False:
            print "Enhanced image does not exists. Enhancing."
            im = Image.open(file['image'])
            im = im.convert("L")  # covert to grayscale
            # gabor filter
            angles = utils.calculate_angles(im, W, f, g)
            print "calculating orientation done"
            angles = utils.smooth_angles(angles)
            print "smoothing angles done"
            im = gabor(im, W, angles)
            # im.show()

            # thinning
            im = make_thin(im)
            # im.show()
            if args.save:
                im.save(image_enhanced_loc, "GIF")
            if args.preprocess:
                continue
        else:
            print "Processing enhanced image."
            im = Image.open(image_enhanced_loc)
            im = im.convert("L")  # covert to grayscale
        # continue
        # get image bounds
        (maxx, maxy) = im.size
        bounds = (0, maxx, 0, maxy)
        print "bounds: ", bounds

        angles = utils.calculate_angles(im, W, f, g)
        if args.smooth:
            angles = utils.smooth_angles(angles)

        if i in ffirst_singularity:
            singularity_type = ffirst_singularity[i]
        else:
            singularity_type = None
        result, cores = calculate_singularities(im, angles, int(args.tolerance[0]), W,
                                                singularity_type=singularity_type)
        print "cores: ", cores
        core = chooseCore(cores, W)
        print "Selected core: ", core
        if j == 1:
            ffirst_singularity[i] = core[3]
        # result.show()

        # align bounds with respect to the core
        # bounds2 = alignBounds(core, bounds)

        # calculate minutias
        im, minutiaes = calculate_minutiaes(im)
        im.show()
        print "minutias: ", minutiaes
        # plotAndSaveMinutiae(minutiaes, image_enhanced_minutiae_loc)
        plotAndSaveMinutiasPIL(minutiaes, image_enhanced_minutiae_loc, im.size)

        minutiaes = alignMinutias(core, minutiaes)  # aligning
        print "aminutia: ", minutiaes
        # plotAndSaveMinutiae(minutiaes, image_enhanced_minutiae_aligned_loc)
        plotAndSaveMinutiasPIL(minutiaes, image_enhanced_minutiae_aligned_loc, im.size)

        des = {'minutiaes': minutiaes, 'core': core, 'bounds': bounds}
        print os.path.splitext(file['image'])[0]
        dumpfilename = "/Users/Admin/fvs/samples/dumps2/" + os.path.basename(file['image']) + ".yaml"
        with open(dumpfilename, 'wb') as handle:
            json.dump(des, handle)
        print ""
        if im is not None:
            im.close()

        # break
        # im.show()

def main2():
    # alignMinutias((2,8, 'loop'), [[40, 120, 'bip']])
    parser = argparse.ArgumentParser(
        description="Gabor, thinning, extract minutias, calculate poincare and align minutias.")
    parser.add_argument("image", nargs=1, help="Path to image")
    # parser.add_argument("block_size", nargs=1, help = "Block size")
    parser.add_argument("tolerance", nargs=1, help="Tolerance for Poincare index")
    parser.add_argument("block_size", nargs=1, help="Blocksize")
    parser.add_argument('--preprocess', "-p", action='store_true', help="Preprocess the image with: Gabor filtering and thinning")
    parser.add_argument('--smooth', "-s", action='store_true', help="Use Gauss for smoothing")
    parser.add_argument("--save", action='store_true', help="Save result image as src_image_enhanced.gif")
    args = parser.parse_args()

    imagepath = args.image[0]
    im = None
    #im = Image.open(args.image[0])
    # im = Image.open("ppf-test_enhanced_thinned.gif")
    #im = im.convert("L")  # covert to grayscale
    # im.show()

    W = int(args.block_size[0])
    #W = 16
    print "Block-size: ", W
    singularity_type = None

    f = lambda x, y: 2 * x * y
    g = lambda x, y: x ** 2 - y ** 2

    #print "Gabor filter and thinning, pointcare, extract minutias, alignment, exporting minutias to file ..."
    print "image: ", os.path.splitext(imagepath)
    base_image_name = os.path.splitext(os.path.basename(imagepath))[0]
    image_enhanced_loc = os.path.splitext(imagepath)[0] + "_enhanced_thinned.gif"
    image_enhanced_minutiae_loc = os.path.splitext(imagepath)[0] + "_enhanced_thinned_minutiaes.png"
    image_enhanced_minutiae_aligned_loc = os.path.splitext(imagepath)[0] + "_enhanced_thinned_minutiaes_aligned.png"
    dumpfilename = os.path.splitext(imagepath)[0] + ".yaml"
    print "image enhanced: ", image_enhanced_loc
    # f['image']
    if args.preprocess or os.path.isfile(image_enhanced_loc) == False:
        print "Enhanced image does not exists. Enhancing."
        print "Gabor filtering and thinning ..."
        im = Image.open(imagepath)
        im = im.convert("L")  # covert to grayscale
        # gabor filter
        angles = utils.calculate_angles(im, W, f, g)
        print "calculating orientation done"
        angles = utils.smooth_angles(angles)
        print "smoothing angles done"
        im = gabor(im, W, angles)
        # im.show()

        # thinning
        im = make_thin(im)
        # im.show()
        if args.save:
            im.save(image_enhanced_loc, "GIF")
        if args.preprocess:
            im.close()
            return
    else:
        print "Processing enhanced image."
        if im is not None:
            im.close()
    im = Image.open(image_enhanced_loc)
    im = im.convert("L")  # covert to grayscale
    # continue
    # get image bounds
    (maxx, maxy) = im.size
    bounds = (0, maxx, 0, maxy)
    print "bounds: ", bounds

    angles = utils.calculate_angles(im, W, f, g)
    if args.smooth:
        angles = utils.smooth_angles(angles)

    singularity_type = None
    result, cores = calculate_singularities(im, angles, int(args.tolerance[0]), W,
                                            singularity_type=singularity_type)
    print "cores: ", cores
    core = chooseCore(cores, W)
    print "Selected core: ", core
    # result.show()

    # align bounds with respect to the core
    # bounds2 = alignBounds(core, bounds)

    # calculate minutias
    im, minutiaes = calculate_minutiaes(im)
    im.show()
    print "minutias: ", minutiaes
    # plotAndSaveMinutiae(minutiaes, image_enhanced_minutiae_loc)
    plotAndSaveMinutiasPIL(minutiaes, image_enhanced_minutiae_loc, im.size)

    minutiaes = alignMinutias(core, minutiaes)  # aligning
    print "aminutia: ", minutiaes
    # plotAndSaveMinutiae(minutiaes, image_enhanced_minutiae_aligned_loc)
    plotAndSaveMinutiasPIL(minutiaes, image_enhanced_minutiae_aligned_loc, im.size)

    des = {'minutiaes': minutiaes, 'core': core, 'bounds': bounds}
    with open(dumpfilename, 'wb') as handle:
        json.dump(des, handle)
    print ""
    # break
    # im.show()

if __name__ == "__main__":
    #main1()
    main2()