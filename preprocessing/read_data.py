import os, struct
from array import array as pyarray
from PIL import Image
import numpy as np


def resize_image(original_image):
    width, height = 16, 16
    resize_image = np.zeros(shape=(width, height))

    for W in range(width):
        for H in range(height):
            new_width = int(W * original_image.shape[0] / width)
            new_height = int(H * original_image.shape[1] / height)
            resize_image[W][H] = original_image[new_width][new_height]

    return resize_image


def load_mnist(dataset="training", digits=np.arange(10), path=".", no_of_imgs=60000):
    if dataset == "training":
        fname_img = os.path.join(path, 'train-images-idx3-ubyte')
        fname_lbl = os.path.join(path, 'train-labels-idx1-ubyte')
    elif dataset == "testing":
        fname_img = os.path.join(path, 't10k-images.idx3-ubyte')
        fname_lbl = os.path.join(path, 't10k-labels.idx1-ubyte')

    else:
        raise ValueError("dataset must be 'testing' or 'training'")

    print("Reading images")
    flbl = open(fname_lbl, 'rb')
    magic_nr, size = struct.unpack(">II", flbl.read(8))
    lbl = pyarray("b", flbl.read())
    flbl.close()

    fimg = open(fname_img, 'rb')
    magic_nr, size, rows, cols = struct.unpack(">IIII", fimg.read(16))
    img = pyarray("B", fimg.read())
    fimg.close()

    ind = [k for k in range(size) if lbl[k] in digits]
    N = size  # int(len(ind) * size/100.)
    images = np.zeros((no_of_imgs, 16, 16), dtype=np.uint8)
    labels = np.zeros((no_of_imgs, 1), dtype=np.int8)
    for i in range(no_of_imgs):  # int(len(ind) * size/100.)):
        temp = np.array(img[ind[i] * rows * cols: (ind[i] + 1) * rows * cols]).reshape((rows, cols))
        images[i] = resize_image(temp)
        labels[i] = lbl[ind[i]]
    labels = [label[0] for label in labels]
    size = images.shape
    images = images.reshape((size[0], size[1] * size[2]))

    return  images / 255 # (images / 255.0) * 0.99 + 0.01


def load_fashion_mnist(dataset="training", digits=np.arange(10), path=".", no_of_imgs=60000):
    if dataset == "training":
        fname_img = os.path.join(path, 'Fashion-mnist\\train-images-idx3-ubyte')  # original path: 'train-images-id5k-ubyte'
        fname_lbl = os.path.join(path, 'Fashion-mnist\\train-labels-idx1-ubyte')  # original path: 'train-labels-id5k-ubyte'
    # elif dataset == "testing":
    # fname_img = os.path.join(path, 't10k-images.idx3-ubyte')
    # fname_lbl = os.path.join(path, 't10k-labels.idx1-ubyte')

    else:
        raise ValueError("dataset must be 'testing' or 'training'")

    print("Reading images")
    flbl = open(fname_lbl, 'rb')
    magic_nr, size = struct.unpack(">II", flbl.read(8))
    lbl = pyarray("b", flbl.read())
    flbl.close()

    fimg = open(fname_img, 'rb')
    magic_nr, size, rows, cols = struct.unpack(">IIII", fimg.read(16))
    img = pyarray("B", fimg.read())
    fimg.close()

    ind = [k for k in range(size) if lbl[k] in digits]
    N = size  # int(len(ind) * size/100.)
    images = np.zeros((no_of_imgs, 16, 16), dtype=np.uint8)
    labels = np.zeros((no_of_imgs, 1), dtype=np.int8)
    for i in range(no_of_imgs):  # int(len(ind) * size/100.)):
        temp = np.array(img[ind[i] * rows * cols: (ind[i] + 1) * rows * cols]).reshape((rows, cols))
        images[i] = resize_image(temp)
        labels[i] = lbl[ind[i]]
    labels = [label[0] for label in labels]
    size = images.shape
    images = images.reshape((size[0], size[1] * size[2]))

    return images / 255


def load_image(path):
    im = Image.open(path).convert('L')
    # im = im.resize((32, 32))
    im = im.resize((16, 16))
    im = np.array(im)
    return im


def load_yale():
    print("Reading Images")
    dir_path = "./Yale"
    # label_dict = {}
    # for i in range(15):
    # 	label_dict[str(i+1)] = i
    image_list = []
    for filename in sorted(os.listdir(dir_path)):
        im = load_image(os.path.join(dir_path, filename))
        image_list.append(im)
    image_list = np.array(image_list)
    size = image_list.shape
    image_list = image_list.reshape((size[0], size[1] * size[2]))

    # return image_list[:50]/5000 
    return image_list / 3000 # 3000


def wine_data(type="red"):
    N = 450
    if type == "red":
        file = open("winequality-red.csv", "r")
    else:
        file = open("winequality-white.csv", "r")
        N = 900
    X = []
    line = file.readline()  # Igonre heading
    line = file.readline()
    while line:
        line = [float(i) for i in line.split(";")]
        X.append(line[:-1])
        line = file.readline()
    X = np.array(X)
    return X / N


def air_quality_data():
    file = open("AirQualityUCI.csv", "r")
    X = []
    line = file.readline()  # Ignore heading
    line = file.readline()

    while line:
        line = line.split(";")[2:-2]
        if line[0] == '':
            break
        temp = []

        for l in line:
            comma = -1
            for i in range(len(l)):
                if l[i] == ',':
                    comma = i
                    break
            if comma == -1:
                temp.append(float(l))
            else:
                temp.append(float(l[:comma]))

        X.append(temp)
        line = file.readline()
    X = np.array(X)
    for i in range(len(X[0])):
        if np.max(abs(X[:, i])) > 200:
            X[:, i] = X[:, i] / 2000

    return X / 3500


def parkinson_data():
    file = open("parkinsons_updrs.data", 'r')
    X = []

    line = file.readline()  # Ignore headers
    line = file.readline()

    while line:
        line = [float(i) for i in line.split(",")]
        line = line[1:8] + line[9:18]
        X.append(line)
        line = file.readline()

    X = np.array(X)

    return X / 1300
