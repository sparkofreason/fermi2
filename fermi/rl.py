import numpy
from scipy.optimize import minimize

def image_dot(i1, i2):
    return sum(sum(numpy.multiply(i1, i2)))

def image_conv(z, psfs):
    r = numpy.zeros(psfs[0].shape)
    for i in range(len(z)):
        r = numpy.add(r, z[i]*psfs[i])
    return r

def rl1(psfs, z, data):
    im0 = image_conv(z, psfs)
    rat = numpy.divide(data, im0)
    zk = [0]*len(z)
    for i in range(len(z)):
        zk[i] = z[i]*image_dot(rat, psfs[i])
    return zk

def loss_map(psfs, z, data):
    im0 = image_conv(z, psfs)
    return numpy.subtract(im0, numpy.multiply(data, numpy.log(im0)))
def d2_loss_map(psfs, z, data):
    im0 = image_conv(z, psfs)
    return numpy.square(numpy.subtract(data, im0))
def loss(psfs, z, data, loss_map):
    return sum(sum(loss_map(psfs, z, data)))
