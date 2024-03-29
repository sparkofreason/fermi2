import itertools
import math
from mmappickle.dict import mmapdict
import multiprocessing
import numpy
from scipy.stats import poisson, skellam, norm


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
        p[small_left] = skellam.cdf(k[small_left], mu1[small_left], mu2[small_left])

    if p[small_right].size > 0:
        ksr = k[small_right]
        mu1sr = mu1[small_right]
        mu2sr = mu2[small_right]
        p1 = skellam.sf(ksr, mu1sr, mu2sr)
        p2 = skellam.pmf(ksr, mu1sr, mu2sr)
        p[small_right] =  p1 + p2 

    if p[large_left].size > 0:
        p[large_left] = norm.cdf(
            k[large_left], diff_mu[large_left], sigma_mu[large_left]
        )

    if p[large_right].size > 0:
        p[large_right] = norm.sf(
            k[large_right],
            diff_mu[large_right],
            sigma_mu[large_right],
        )

    return numpy.reshape(p, k.shape)


def poisson_tail(x, mu):
    p = numpy.ones(x.shape)
    left = x <= mu
    right = x > mu
    p[left] = poisson.cdf(x[left], mu[left])
    p[right] = poisson.pmf(x[right], mu[right]) + poisson.sf(x[right], mu[right])
    return numpy.reshape(p, x.shape)


def sidak(alpha, j):
    return 1 - (1 - alpha) ** (1 / (2 ** (2 * j)))


def pvalues1(*args):
    if (len(args) == 6):
        k, k_model, mu1, mu2, j, direction = args
        print(direction, j, "******************************")
        return skellam_tail(k, mu1, mu2)
    else:
        a_counts, a_model = args
        return poisson_tail(a_counts, a_model)


def threshold1_impl(k, k_model, p, alpha_j, j, direction):
    print(direction, j, "******************************")
    print("alpha_j", alpha_j)
    mask = p < alpha_j/2
    k[mask] = k[mask] - k_model[mask]
    k[~mask] = 0
    return k


def threshold1(*args):
    if (len(args) == 6):
        return threshold1_impl(*args)
    else:
        a_counts, a_model, ap, alpha_j = args
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


def pvalues_pool(args, poolWorkers = None):
    phs = []
    pvs = []
    pds = []

    if (poolWorkers == None):
        # TODO - This is broken.
        ks = map(pvalues1, args)
    else:
        ctx = multiprocessing.get_context("fork")
        p = ctx.Pool(poolWorkers, maxtasksperchild=1)
        with p:
            pvals = []
            res = []
            for a in args:
                res.append(p.apply_async(pvalues1, a))
            for r in res:
                pvals.append(r.get())
    
    for pv in chunked_iterable(pvals, 3):
        if (len(pv) == 3):
            phs.append(pv[0])
            pvs.append(pv[1])
            pds.append(pv[2])
        else:
            pa = pv[0]

    return {'pvalue': {'a': pa, 'hs': phs, 'vs': pvs, 'ds': pds}}


def threshold(args):
    hs = []
    vs = []
    ds = []

    ks = []
    for a in args:
        ks.append(threshold1(*a))
    
    for k in chunked_iterable(ks, 3):
        if (len(k) == 3):
            hs.append(k[0])
            vs.append(k[1])
            ds.append(k[2])
        else:
            a = k[0]

    return {'wavelet': {'a': a, 'hs': hs, 'vs': vs, 'ds': ds}}


# def haar_threshold(counts, model, alpha, jmin=0, fwer=None):
#     args = skellam_inputs(counts, model, alpha, jmin, fwer)
#     return haar_threshold_pool(args)


# def haar_threshold_sphere(counts, model, alpha, jmin=0, fwer=None):
#     a_counts = spherify(counts)
#     a_model = spherify(model)
#     return haar_threshold(a_counts, a_model, alpha, jmin, fwer)


def inv_haar_j(J, j, a, h, v, d):
    s = 2 ** (J - j - 1)
    ah = a + h 
    av = a + v 
    ad = a + d 
    vd = v + d 
    hd = h + d 
    hv = h + v 
    xx = ah + vd 
    xr = av - hd 
    dx = ah - vd 
    dr = ad - hv 
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
def inv_haar_level(J, j, a, h, v, d, direction=None ):
    j = j - 1
    zs = 0 #numpy.zeros(a.shape)
    hsp = list(itertools.repeat(zs, J))
    vsp = list(itertools.repeat(zs, J))
    dsp = list(itertools.repeat(zs, J))
    if (j == J):
        ap = a
    else:
        ap = numpy.zeros(a.shape)
        if direction == None:
            hsp[j] = h
            vsp[j] = v
            dsp[j] = d
        elif direction == 'hs':
            hsp[j] = h
        elif direction == 'vs':
            vsp[j] = v
        elif direction == 'ds':
            dsp[j] = d
        
    return inv_haar_sphere({'wavelet': {'a': ap, 'hs': hsp, 'vs': vsp, 'ds': dsp}})

def level_key(*keys):
    return '/'.join(map(str, keys))


def pvalue_inputs(counts, model, jmin=0,):
    a_counts = numpy.copy(counts)
    a_model = numpy.copy(model)

    rows, cols = a_counts.shape

    J = int(math.ceil(math.log(max(rows, cols), 2)))

    for j in range(J - 1, jmin - 1, -1):
        a_counts, h_counts, v_counts, d_counts = haar_1(haar_sums_j(J, j, a_counts))

        sums_model = haar_sums_j(J, j, a_model)
        a_model, h_model, v_model, d_model = haar_1(sums_model)

        (h1, h2), (v1, v2), (d1, d2) = sums_model

        yield (h_counts, h_model, h1, h2, j, "Horizontal")   
        yield (v_counts, v_model, v1, v2, j, "Vertical")
        yield (d_counts, d_model, d1, d2, j, "Diagonal")
    yield (a_counts, a_model)


def run_pvalues(count_data, total_model, jmin, filename, poolWorkers = 1):
    args = pvalue_inputs(spherify(count_data), spherify(total_model), jmin)
    
    result = pvalues_pool(args, poolWorkers=poolWorkers)

    m = mmapdict(filename)
    for j in range(1, 13):
        for d in ['hs', 'vs', 'ds']:
            m[level_key('pvalue', d, j)] = result['pvalue'][d][j-1]
    m[level_key('pvalue', 'a')] = result['pvalue']['a']
    m['count_data'] = count_data
    m['total_model'] = total_model
    m.vacuum()


def threshold_inputs(pvs, alpha_fn, jmin=0):
    a_counts = spherify(pvs['count_data'])
    a_model = spherify(pvs['total_model'])

    rows, cols = a_counts.shape

    J = int(math.ceil(math.log(max(rows, cols), 2)))

    for j in range(J - 1, jmin - 1, -1):
        alpha_j = alpha_fn(j)

        a_counts, h_counts, v_counts, d_counts = haar_1(haar_sums_j(J, j, a_counts))

        a_model, h_model, v_model, d_model = haar_1(haar_sums_j(J, j, a_model))

        yield (h_counts, h_model, pvs[level_key('pvalue', 'hs', J-j)], alpha_j, j, "Horizontal")   
        yield (v_counts, v_model, pvs[level_key('pvalue', 'vs', J-j)], alpha_j, j, "Vertical")
        yield (d_counts, d_model, pvs[level_key('pvalue', 'ds', J-j)], alpha_j, j, "Diagonal")
    yield (a_counts, a_model, pvs[level_key('pvalue', 'a')], alpha_j)


def run_threshold(pvalues_filename, alpha_fn, jmin, filename, level_recs = False, poolWorkers = 1):
    pvs = mmapdict(pvalues_filename, readonly=True)
    args = threshold_inputs(pvs, alpha_fn, jmin)
    
    result = threshold(args)

    m = mmapdict(filename)
    count_rec = inv_haar_sphere(result)
    m['count_rec'] = count_rec
    if level_recs:
        ctx = multiprocessing.get_context("fork")
        p = ctx.Pool(poolWorkers, maxtasksperchild=1)
        with p:
            res = []
            w = result['wavelet']
            J = len(w['hs'])
            for j in range(1, J + 1):
                for d in ['hs', 'vs', 'ds']:
                    m[level_key('wavelet', d, j)] = result['wavelet'][d][j-1]
            for j in range(1, J + 2):
                print('Reconstructing level: ', j)
                jj = j-1 if j < J else J - 1
                res.append(p.apply_async(inv_haar_level, (J, j, w['a'], w['hs'][jj], w['vs'][jj], w['ds'][jj])))
            for j in range(1, J + 2):
                rec = res[j-1].get()
                m[level_key('level_recs', j)] = rec
            m[level_key('wavelet', 'a')] = result['wavelet']['a']

    m.vacuum()
