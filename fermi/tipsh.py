import itertools
import math
import multiprocessing
import numpy
from scipy.stats import poisson, skellam
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

numpy2ri.activate()

stats = importr("stats")
skel = importr("skellam")

# Really just a 2D numpy.roll
def roll_sphere(arr, lat_roll, lon_roll):
    rolled = numpy.roll(arr, lon_roll, 1)
    rolled = numpy.roll(rolled, lat_roll, 0)
    return rolled


# Given arr that represents a cylindrical projection of a sphere,
# stack a shifted and flipped version such that going over either
# effectively takes you to the "other side". So, for example, if 
# you were following longitude 30 degrees, and went up through
# latitude 89->90->89, you would land at longitude 30+180 = 210.
def spherify(arr):
    _, cols = arr.shape
    spharr = numpy.vstack((arr, numpy.flipud(numpy.roll(arr, cols // 2, 1))))
    return spharr


# Undo spherify and return the original array
def unspherify(spharr):
    rows, _ = spharr.shape
    return spharr[0 : rows // 2, :]


# Same as unspherify, but operating on the to stacked array.
# Handy for debugging
def unspherify_top(spharr):
    rows, cols = spharr.shape
    return numpy.flipud(numpy.roll(spharr[rows // 2 : rows, :], -cols // 2, 1))


# Returns shifted versions of arr suitable for computing the Haar
# transform at a given scale.
# Return:
# - xx - original arr
# - xr 
def haar_vals(J, j, arr):
    s = 2 ** (J - j - 1)
    xx = arr
    xr = roll_sphere(arr, 0, -s)
    dx = roll_sphere(arr, -s, 0)
    dr = roll_sphere(arr, -s, -s)
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
    large = (mu1 >= 100) | (mu2 > 100)
    small = ~large
    left = k <= diff_mu
    right = ~left
    small_left = small & left
    small_right = small & right
    large_left = large & left
    large_right = large & right

    if p[small_left].size > 0:
        print("small_left")
        p[small_left] = skellam.cdf(k[small_left], mu1[small_left], mu2[small_left])
        # p[small_left] = skel.pskellam(
        #     k[small_left], mu1[small_left], mu2[small_left], lower_tail=True
        # )

    if p[small_right].size > 0:
        print("small_right")
        ksr = k[small_right]
        mu1sr = mu1[small_right]
        mu2sr = mu2[small_right]
        # try:
        #     p2 = skel.dskellam(ksr, mu1sr, mu2sr)
        # except:
        #     p2 = skel.dskellam_sp(ksr, mu1sr, mu2sr)
        # p1 = skel.pskellam(ksr, mu1sr, mu2sr, lower_tail=False)
        p1 = skellam.sf(ksr, mu1sr, mu2sr)
        p2 = skellam.pmf(ksr, mu1sr, mu2sr)
        p[small_right] =  p1 + p2 

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

    return numpy.reshape(p, k.shape)


def poisson_tail(x, mu):
    p = numpy.ones(x.shape)
    left = x <= mu
    right = x > mu
    p[left] = stats.ppois(x[left], mu[left])
    p[right] = stats.ppois(x[right], mu[right], lower_tail=False)
    return numpy.reshape(p, x.shape)

def sidak(alpha, j):
    return 1 - (1 - alpha) ** (1 / (2 ** (2 * j)))

def skellam_inputs(counts, model, alpha, jmin=0, fwer=None):
    a_counts = numpy.copy(counts)
    a_model = numpy.copy(model)

    rows, cols = a_counts.shape

    J = int(math.ceil(math.log(max(rows, cols), 2)))

    for j in range(J - 1, jmin - 1, -1):
        f = 2 ** (J - j)
        if fwer == "sidak":
            alpha_j = sidak(alpha, j)
        elif fwer == "bonferroni":
            alpha_j = alpha / 2 ** (2 * j)
        else:
            alpha_j = alpha

        print("sums")
        a_counts, h_counts, v_counts, d_counts = haar_1(haar_sums_j(J, j, a_counts))

        sums_model = haar_sums_j(J, j, a_model)
        a_model, h_model, v_model, d_model = haar_1(sums_model)

        (h1, h2), (v1, v2), (d1, d2) = sums_model
        yield (h_counts, h_model, h1, h2, alpha_j, j, "Horizontal")   
        yield (v_counts, v_model, v1, v2, alpha_j, j, "Vertical")
        yield (d_counts, d_model, d1, d2, alpha_j, j, "Diagonal")
    yield (a_counts, a_model, alpha_j)


def threshold1_impl(k, k_model, mu1, mu2, alpha_j, j, direction):
    print(direction, j, "******************************")
    print("alpha_j", alpha_j)
    p = skellam_tail(k, mu1, mu2)
    mask = p < alpha_j/2
    k[mask] = k[mask] - k_model[mask]
    k[~mask] = 0
    return k, p


def threshold1(*args):
    if (len(args) == 7):
        return threshold1_impl(*args)
    else:
        a_counts, a_model, alpha_j = args
        ap = poisson_tail(a_counts, a_model)
        a_mask = ap < alpha_j / 2
        a_counts[a_mask] = a_counts[a_mask] - a_model[a_mask]
        a_counts[~a_mask] = 0
        return a_counts, ap


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


def haar_threshold_pool(args, poolWorkers = None):
    hs = []
    vs = []
    ds = []
    phs = []
    pvs = []
    pds = []

    if (poolWorkers == None):
        # TODO - This is broken.
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
            hs.append(k[0][0])
            phs.append(k[0][1])
            vs.append(k[1][0])
            pvs.append(k[1][1])
            ds.append(k[2][0])
            pds.append(k[2][1])
        else:
            a = k[0][0]
            pa = k[0][1]

    return {'wavelet': {'a': a, 'hs': hs, 'vs': vs, 'ds': ds},
            'pvalue': {'a': pa, 'hs': phs, 'vs': pvs, 'ds': pds}}


def haar_threshold(counts, model, alpha, jmin=0, fwer=None):
    args = skellam_inputs(counts, model, alpha, jmin, fwer)
    return haar_threshold_pool(args)


def haar_threshold_sphere(counts, model, alpha, jmin=0, fwer=None):
    a_counts = spherify(counts)
    a_model = spherify(model)
    return haar_threshold(a_counts, a_model, alpha, jmin, fwer)


def inv_haar_j(J, j, a, h, v, d):
    s = 2 ** (J - j - 1)
    ah = a + h #numpy.add(a, h)
    av = a + v #numpy.add(a, v)
    ad = a + d #numpy.add(a, d)
    vd = v + d #numpy.add(v, d)
    hd = h + d #numpy.add(h, d)
    hv = h + v #numpy.add(h, v)
    xx = ah + vd #numpy.add(ah, vd)
    xr = av - hd #numpy.subtract(av, hd)
    dx = ah - vd #numpy.subtract(ah, vd)
    dr = ad - hv #numpy.subtract(ad, hv)
    arr = (
        xx + roll_sphere(xr, 0, s) + roll_sphere(dx, s, 0) + roll_sphere(dr, s, s)
    ) / 16 # 4 terms for each shift times 4 values in summed in each coarse coefficient
    return arr


def inv_haar(result):
    wavelet = result['wavelet']
    a = wavelet['a']
    hs = wavelet['hs']
    vs = wavelet['vs']
    ds = wavelet['ds']
    rows, cols = a.shape
    J = int(math.ceil(math.log(max(rows, cols), 2)))
    nj = len(hs)
    for i in range(nj - 1, -1, -1):
        j = J - i - 1
        a = inv_haar_j(J, j, a, hs[i], vs[i], ds[i])
    return a


def inv_haar_sphere(result):
    spharr = inv_haar(result)
    return (unspherify(spharr) + unspherify_top(spharr))/2

# I'm missing something here, seems like this should be a lot faster
def inv_haar_level(j, result, direction=None ):
    w = result['wavelet']
    j = j - 1
    zs = numpy.zeros(w['a'].shape)
    hsp = list(itertools.repeat(zs, len(w['hs'])))
    vsp = list(itertools.repeat(zs, len(w['vs'])))
    dsp = list(itertools.repeat(zs, len(w['ds'])))
    ap = zs
    if direction == None:
        hsp[j] = w['hs'][j]
        vsp[j] = w['vs'][j]
        dsp[j] = w['ds'][j]
    elif direction == 'hs':
        hsp[j] = w['hs'][j]
    elif direction == 'vs':
        vsp[j] = w['vs'][j]
    elif direction == 'ds':
        dsp[j] = w['ds'][j]
        
    return inv_haar_sphere({'wavelet': {'a': ap, 'hs': hsp, 'vs': vsp, 'ds': dsp}})