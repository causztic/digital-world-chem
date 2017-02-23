import numpy
from laguerre import normalized_radial_solution
from legendre import normalized_angular_solution

def hydrogen_wave_func(n, l, m, roa, nx, ny, nz):
    xx, yy, zz = numpy.meshgrid(nx, ny, nz)
    normalized_angular_solution(m,l)

hydrogen_wave_func(2,1,1,8,2,2,2)
