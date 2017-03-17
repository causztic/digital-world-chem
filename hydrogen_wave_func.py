import math
import numpy
import scipy.constants as c
from functions import cartesian_to_spherical
from laguerre import normalized_radial_solution, radial_wave_func
from legendre import normalized_angular_solution, angular_wave_func

def hydrogen_wave_func(n, l, m, roa, nx, ny, nz):
    roa = float(roa)
    xx, yy, zz = numpy.meshgrid(numpy.linspace(-roa, roa, nx),
                                numpy.linspace(-roa, roa, ny), numpy.linspace(-roa, roa, nz))

    vector_r = numpy.vectorize(cartesian_to_spherical)
    v_angularwf = numpy.vectorize(angular_wave_func)
    v_radialwf = numpy.vectorize(radial_wave_func)

    bohr_radius = c.physical_constants['Bohr radius'][0]
    r, theta, phi = vector_r(xx, yy, zz)
    r = r * bohr_radius

    radial_part = v_radialwf(n, l, r)
    angular_part = v_angularwf(m, l, theta, phi)

    rwf = (radial_part * angular_part) ** 2
    mag = numpy.absolute(rwf)

    return numpy.round(xx,5), numpy.round(yy,5), numpy.round(zz,5), numpy.round(mag,5)