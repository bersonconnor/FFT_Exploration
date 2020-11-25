import math
import random


# Computes the real part of the exponential function
#
# @param exp - the exponent of the exponential function
def get_real_exp(exp):
    return math.cos(exp)


# Computes the imaginary part of the exponential function
#
# @param exp - the exponent of the exponential function
def get_img_exp(exp):
    return math.sin(exp)


# Computes the value of the exponential function
#
# @param exp - the exponent of the exponential function
def CEXPO(exp):
    real = get_real_exp(exp)
    img = get_img_exp(exp)
    return real, img


# Returns the result of multiplying complex numbers
#
# @param z_1 - complex number 1
# @param z_2 - complex number 2
def CMUL(z_1, z_2):
    (x_1, y_1) = z_1
    (x_2, y_2) = z_2
    real = x_1 * x_2 - y_1 * y_2
    img = x_1 * y_2 + x_2 * y_1
    result = (real, img)
    return result


# Returns the result of adding complex numbers
#
# @param z_1 - complex number 1
# @param z_2 - complex number 2
def CADD(z_1, z_2):
    (x_1, y_1) = z_1
    (x_2, y_2) = z_2
    real = x_1 + x_2
    img = y_1 + y_2
    result = (real, img)
    return result

# Returns the result of subtracting complex numbers
#
# @param z_1 - complex number 1
# @param z_2 - complex number 2
def CSUB(z_1, z_2):
    (x_1, y_1) = z_1
    (x_2, y_2) = z_2
    real = x_1 - x_2
    img = y_1 - y_2
    result = (real, img)
    return result


# Round the real and imaginary parts of a complex number
# to the closest integer
#
# @param z - the complex number we are rounding
def round_img(z):
    (x, y) = z
    (u, v) = (int(round(x)), int(round(y)))
    z = (u, v)
    return z


# Returns a randomly generated imaginary number
# that is no greater than the given maximum
#
# @param max_val - largest real value
def gen_img(max_val):
    real = random.randint(0, max_val)
    img = random.randint(0, max_val)
    z = (real, img)
    return z


# Returns a randomly generated list imaginary numbers
# that is length size with no number larger than max_val
#
# @param max_val - largest real value
# @param power - length of returned list as 2^power
def gen_img_list(max_val, power):
    res = []
    for i in range(2 ** power):
        res.append(gen_img(max_val))
    return res
