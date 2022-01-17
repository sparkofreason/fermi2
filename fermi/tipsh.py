import gc
import itertools
import math
import multiprocessing
import numpy
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

numpy2ri.activate()

stats = importr("stats")
skel = importr("skellam")

def roll_sphere(arr, lat_roll, lon_roll):
    rolled = numpy.roll(arr, lon_roll, 1)
    rolled = numpy.roll(rolled, lat_roll, 0)
    return rolled

def spherify(arr):
    _, cols = arr.shape
    spharr = numpy.vstack((arr, numpy.flipud(numpy.roll(arr, cols // 2, 1))))
    return spharr


def unspherify(spharr):
    rows, _ = spharr.shape
    return spharr[0 : rows // 2, :]


def haar_vals(J, j, arr):
    s = 2 ** (J - j - 1)
    xx = arr
    xr = roll_sphere(arr, 0, s)
    dx = roll_sphere(arr, s, 0)
    dr = roll_sphere(arr, s, s)
    return (xx, xr, dx, dr)


def haar_sums_j(J, j, arr):
    xx, xr, dx, dr = haar_vals(J, j, arr)
    h1 = numpy.add(xx, dx)
    h2 = numpy.add(xr, dr)
    v1 = numpy.add(xx, xr)
    v2 = numpy.add(dx, dr)
    d1 = numpy.add(xx, dr)
    d2 = numpy.add(xr, dx)
    return ((h1, h2), (v1, v2), (d1, d2))


def haar_1(sums):
    (h1, h2), (v1, v2), (d1, d2) = sums
    a = numpy.add(h1, h2) 
    h = numpy.subtract(h1, h2) 
    v = numpy.subtract(v1, v2) 
    d = numpy.subtract(d1, d2) 
    return (a, h, v, d)


def haar(a, jmin=0):
    rows, cols = a.shape
    J = int(math.ceil(math.log(max(rows, cols), 2)))
    hs = []
    vs = []
    ds = []
    for j in range(J - 1, jmin - 1, -1):
        sums = haar_sums_j(J, j, a)
        a, h, v, d = haar_1(sums)
        hs.append(h)
        vs.append(v)
        ds.append(d)
    return a, hs, vs, ds


def haar_sphere(arr, jmin=0):
    a = spherify(arr)
    return haar(a, jmin)


def skellam_tail(k, mu1, mu2):
    p = numpy.ones(k.shape)
    diff_mu = numpy.subtract(mu1, mu2)
    sigma_mu = numpy.sqrt(numpy.add(mu1, mu2))
    large = (mu1 >= 1000) | (mu2 > 1000)
    small = ~large
    left = k <= diff_mu
    right = ~left
    small_left = small & left
    small_right = small & right
    large_left = large & left
    large_right = large & right

    if p[small_left].size > 0:
        print("small_left")
        p[small_left] = skel.pskellam(
            k[small_left], mu1[small_left], mu2[small_left], lower_tail=True
        )

    if p[small_right].size > 0:
        print("small_right")
        p[small_right] = skel.pskellam(
            k[small_right], mu1[small_right], mu2[small_right], lower_tail=False
        )

    if p[large_left].size > 0:
        print("large_left")
        p[large_left] = stats.pnorm(
            k[large_left], diff_mu[large_left], sigma_mu[large_left], lower_tail=True
        )

    if p[large_right].size > 0:
        print("large_right")
        p[large_right] = stats.pnorm(
            k[large_right],
            diff_mu[large_right],
            sigma_mu[large_right],
            lower_tail=False,
        )

    return p


def poisson_tail(x, mu):
    p = numpy.ones(x.shape)
    left = x <= mu
    right = x > mu
    p[left] = stats.ppois(x[left], mu[left])
    p[right] = stats.ppois(x[right], mu[right], lower_tail=False)
    return p

def skellam_inputs(counts, model, alpha, jmin=0, fwer=None, mode=None):
    a_counts = numpy.copy(counts)
    a_model = numpy.copy(model)

    rows, cols = a_counts.shape

    J = int(math.ceil(math.log(max(rows, cols), 2)))

    for j in range(J - 1, jmin - 1, -1):
        f = 2 ** (J - j)
        if fwer == "sidak":
            alpha_j = 1 - (1 - alpha) ** (1 / (2 ** (2 * j)))
        elif fwer == "bonferroni":
            alpha_j = alpha / 2 ** (2 * j)
        else:
            alpha_j = alpha

        print("sums")
        a_counts, h_counts, v_counts, d_counts = haar_1(haar_sums_j(J, j, a_counts))

        sums_model = haar_sums_j(J, j, a_model)
        a_model, h_model, v_model, d_model = haar_1(sums_model)

        (h1, h2), (v1, v2), (d1, d2) = sums_model
        yield (h_counts, h_model, h1, h2, alpha_j, j, "Horizontal", mode)
        yield (v_counts, v_model, v1, v2, alpha_j, j, "Vertical", mode)
        yield (d_counts, d_model, d1, d2, alpha_j, j, "Diagonal", mode)
    yield (a_counts, a_model, alpha_j, mode)

def threshold1_impl(k, k_model, mu1, mu2, alpha_j, j, direction, mode):
    print(direction, j, "******************************")
    print("alpha_j", alpha_j)
    p = skellam_tail(k, mu1, mu2)
    if mode == 'pvalue':
        return p
    mask = p < alpha_j/2
    k[mask] = k[mask] - k_model[mask]
    k[~mask] = 0
    return k

def threshold1(*args):
    if (len(args) == 8):
        return threshold1_impl(*args)
    else:
        a_counts, a_model, alpha_j, mode = args
        ap = poisson_tail(a_counts, a_model)
        if mode == 'pvalue':
            return ap
        a_mask = ap < alpha_j / 2
        a_counts[a_mask] = a_counts[a_mask] - a_model[a_mask]
        a_counts[~a_mask] = 0
        return a_counts

def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk

def haar_threshold_foo(args, poolWorkers = None):
    hs = []
    vs = []
    ds = []

    if (poolWorkers == None):
        ks = map(threshold1, args)
    else:
        ctx = multiprocessing.get_context("fork")
        p = ctx.Pool(poolWorkers, maxtasksperchild=1)
        with p:
            #ks = p.imap(threshold1, args, 1)
            ks = []
            res = []
            for a in args:
            #    _, k, *rest = a
            #    ks.append(k)
                res.append(p.apply_async(threshold1, a))
            for r in res:
                ks.append(r.get())
    
    for k in chunked_iterable(ks, 3):
        if (len(k) == 3):
            hs.append(k[0])
            vs.append(k[1])
            ds.append(k[2])
        else:
            a = k[0]

    return a, hs, vs, ds

def haar_threshold(counts, model, alpha, jmin=0, fwer=None):
    args = skellam_inputs(counts, model, alpha, jmin, fwer)
    return haar_threshold_foo(args)

def haar_threshold_sphere(counts, model, alpha, jmin=0, fwer=None):
    a_counts = spherify(counts)
    a_model = spherify(model)
    return haar_threshold(a_counts, a_model, alpha, jmin, fwer)


def inv_haar_j(J, j, a, h, v, d):
    s = 2 ** (J - j - 1)
    ah = numpy.add(a, h)
    av = numpy.add(a, v)
    ad = numpy.add(a, d)
    vd = numpy.add(v, d)
    hd = numpy.add(h, d)
    hv = numpy.add(h, v)
    xx = numpy.add(ah, vd)
    xr = numpy.subtract(av, hd)
    dx = numpy.subtract(ah, vd)
    dr = numpy.subtract(ad, hv)
    spharr = (
        xx + roll_sphere(xr, 0, -s) + roll_sphere(dx, -s, 0) + roll_sphere(dr, -s, -s)
    ) / 16
    return spharr


def inv_haar(a, hs, vs, ds):
    rows, cols = a.shape
    J = int(math.ceil(math.log(max(rows, cols), 2)))
    nj = len(hs)
    for i in range(nj - 1, -1, -1):
        j = J - i - 1
        a = inv_haar_j(J, j, a, hs[i], vs[i], ds[i])
    return a


def inv_haar_sphere(a, hs, vs, ds):
    return unspherify(inv_haar(a, hs, vs, ds))

# I'm missing something here, seems like this should be a lot faster
def inv_haar_level(j, a, hs, vs, ds):
    j = j - 1
    zs = numpy.zeros(a.shape)
    hsp = list(itertools.repeat(zs, len(hs)))
    vsp = list(itertools.repeat(zs, len(vs)))
    dsp = list(itertools.repeat(zs, len(ds)))
    ap = zs
    hsp[j] = hs[j]
    vsp[j] = vs[j]
    dsp[j] = ds[j]
    return inv_haar_sphere(ap, hsp, vsp, dsp)