import math
import time
from cmath import exp, pi
import img_utils
from img_utils import CEXPO, CMUL, CADD, CSUB


# Computes the value of seq[n] * e^(-i2pi/N * kn)
#
# @param seq - input for DFT
# @param n - current index of sum
# @param k - current index of input
def compute_formula(seq, n, k):
    x_n = seq[n]
    exp = (-2 * pi * n * k) / len(seq)
    exp_tpl = CEXPO(exp)
    result = CMUL(x_n, exp_tpl)

    return result


# Computes the value of summing the formula seq[n] * e^(-i2pi/N * kn)
# while iterating over n with a fixed k
#
# # @param seq - input for DFT
# # @param n - current index of sum
def compute_sum(seq, k):
    result = (0, 0)
    n = 0
    while n < len(seq):
        val = compute_formula(seq, n, k)
        result = CADD(result, val)
        n += 1

    return result


# Returns the result of applying the Discrete Fourier Transform
#
# @param samples - a finite sequence of equally-spaced samples of a function
# @returns a same-length sequence of equally-spaced samples of
# the discrete fourier transform (a complex valued function
# of frequency)
def DFT(samples):
    N = len(samples)
    result = []
    k = 0
    while k < N:
        total = (0, 0)
        n = 0
        while n < N:
            total = CADD(total, CMUL(samples[n], CEXPO(-2 * pi * n * k / N)))
            n += 1

        result.append(total)
        k += 1

    return result


# Runs the recursive Cooley-Turkey Algorithm
#
# @param samples - the input sequence
def CTFFT(samples):
    N = len(samples)
    if N <= 1:
        return samples
    else:
        evens = CTFFT(samples[0::2])
        odds = CTFFT(samples[1::2])

        odds = [CMUL(CEXPO(-2 * pi * k / N), odds[k]) for k in range(N // 2)]
        left = [CADD(evens[k], odds[k]) for k in range(N / 2)]
        right = [CSUB(evens[k], odds[k]) for k in range(N / 2)]
        return left + right


# constant used repeatedly for computation
exp_const = -2 * pi
# Runs the recursive Cooley-Turkey Algorithm but
# cuts off to vanilla DFT as soon as the length
# of the list is less than or equal to cutoff
# with the default being 32
#
# @param samples - the input sequence
def OPT_CTFFT(samples):
    N = len(samples)
    if N <= 1:
        return samples
    else:
        evens = OPT_CTFFT(samples[0::2])
        odds = OPT_CTFFT(samples[1::2])
        odds = [CMUL(CEXPO(exp_const * k / N), odds[k]) for k in range(N // 2)]

        left = []
        right = []
        for k in range(N / 2):
            x_1, y_1 = evens[k]
            x_2, y_2 = odds[k]
            left.append((x_1 + x_2, y_1 + y_2))
            right.append((x_1 - x_2, y_1 - y_2))

        return left + right


# Returns a rounded list of tuples
#
# @param res - list of tuples to round
def round_result(res):
    for i in range(len(res)):
        res[i] = img_utils.round_img(res[i])
    return res


# Runs a random example with both the vanilla DFT and
# the Cooley-Tureky algorithm and times both of their
# speeds
#
# @param max_val - largest real component for samples
# @param power - length of example as 2^power
# @param isCT - true to test vanilla CT and false to not
# @param isCCT - true to test custom CT and false to not
# @param isDFT - true to test DFT and false to not
def test_random_example(max_val, power, isCT=True, isCCT=True, isDFT=True):
    # time variables
    ct_time = 0
    cct_time = 0
    dft_time = 0

    # Generate an example
    example = img_utils.gen_img_list(max_val, power)
    #print("Example:\t\t" + str(example) + "\n")

    if isCT:
        # Run and time the Vanilla Cooley-Turkey algorithm
        ct_start = time.time()
        ct_res = CTFFT(example)
        ct_end = time.time()
        ct_time = ct_end - ct_start
        #print("Cooley-Turkey:\t" + str(ct_time) + " seconds")
        ct_res = round_result(ct_res)
        #print("Cooley-Turkey:\t" + str(ct_res) + "\n")

    if isCCT:
        # Run and time the Cut off Cooley-Turkey algorithm
        cct_start = time.time()
        cct_res = OPT_CTFFT(example)
        cct_end = time.time()
        cct_time = cct_end - cct_start
        #print("Custom CT:\t\t" + str(cct_time) + " seconds")
        cct_res = round_result(cct_res)
        #print("Custom CT:\t\t" + str(cct_res) + "\n")

    if isDFT:
        # Run and time the Vanilla DFT algorithm
        dft_start = time.time()
        dft_res = DFT(example)
        dft_end = time.time()
        dft_time = dft_end - dft_start
        #print("Vanilla DFT:\t" + str(dft_time) + " seconds")
        dft_res = round_result(dft_res)
        #print("Vanilla DFT:\t" + str(dft_res))

    return ct_time, cct_time, dft_time


# Runs num_runs examples on the vanilla CTFFT, custom CTFFT, and
# vanilla DFT algorithms and prints the total time necessary and
# the average time for each
#
# @param max_val - largest real component for samples
# @param power - length of example as 2^power
# @param num_runs - number of examples that will be computed
# @param isCT - true to test vanilla CT and false to not
# @param isCCT - true to test custom CT and false to not
# @param isDFT - true to test DFT and false to not
def test_many_random_examples(max_val, power, num_runs, isCT=True, isCCT=True, isDFT=True):
    # time variables
    ct_total = 0
    cct_total = 0
    dft_total = 0

    worst_ct = -1000000000000000
    worst_cct = -1000000000000000
    worst_dft = -1000000000000000

    best_ct = 1000000000000000
    best_cct = 1000000000000000
    best_dft = 1000000000000000

    i = 0
    while i < num_runs:
        ct_time, cct_time, dft_time = test_random_example(max_val, power, isCT, isCCT, isDFT)
        ct_total += ct_time
        cct_total += cct_time
        dft_total += dft_time

        if ct_time > worst_ct:
            worst_ct = ct_time

        if cct_time > worst_cct:
            worst_cct = cct_time

        if dft_time > worst_dft:
            worst_dft = dft_time

        if ct_time < best_ct:
            best_ct = ct_time

        if cct_time < best_cct:
            best_cct = cct_time

        if dft_time < best_dft:
            best_dft = dft_time

        i += 1

    ct_avg = ct_total / num_runs
    cct_avg = cct_total / num_runs
    dft_avg = dft_total / num_runs
    print("----- Max Value = " + str(max_val) + "; N = 2^" + str(power) + "; Number Runs: " + str(num_runs) + " -----")
    print("Vanilla CT: \n\tTotal Time:\t\t" + str(ct_total) + "\n\tAverage Time:\t" + str(ct_avg) + "\n\tBest Time:\t\t" + str(best_ct) + "\n\tWorst Time:\t\t" + str(worst_ct) + "\n")
    print("Custom CT: \n\tTotal Time:\t\t" + str(cct_total) + "\n\tAverage Time:\t" + str(cct_avg) + "\n\tBest Time:\t\t" + str(best_cct) + "\n\tWorst Time:\t\t" + str(worst_cct) +"\n")
    print("Vanilla DFT: \n\tTotal Time:\t\t" + str(dft_total) + "\n\tAverage Time:\t" + str(dft_avg) + "\n\tBest Time:\t\t" + str(best_dft) + "\n\tWorst Time:\t\t" + str(worst_dft) + "\n\n")


# Runs test_many_random_examples but varies the value of N by
# incrementing the value of the power starting with power = 1
#
# @param max_val - largest real component for samples
# @param min_power - minimum length of examples as 2^power
# @param max_power - maximum length of examples as 2^power
# @param num_runs - number of examples that will be computed
# @param isCT - true to test vanilla CT and false to not
# @param isCCT - true to test custom CT and false to not
# @param isDFT - true to test DFT and false to not
def test_many_examples_with_varied_N(max_val, min_power, max_power, num_runs, isCT=True, isCCT=True, isDFT=True):
    i = min_power
    while i < max_power + 1:
        test_many_random_examples(max_val, i, num_runs, isCT, isCCT, isDFT)
        i += 1


test_many_examples_with_varied_N(10, 1, 12, 200, True, True, False)