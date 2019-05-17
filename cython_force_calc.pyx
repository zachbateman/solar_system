'''
Pulling out expensive calculations to CYTHON
'''
import cython
from functools import lru_cache

# For decorator usage see:  https://stackoverflow.com/questions/4709285/cython-float-division-pyexc-zerodivisionerror-checking?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
@cython.cdivision(True)
# @lru_cache(maxsize=1024)
cpdef calc_composite_force_vector(double grav_const, self, list spacebodies_data):
    cdef double force_vec_x, force_vec_y, force_vec_z
    cdef double composite_x, composite_y, composite_z
    cdef str selfid
    cdef double grav_const_times_mass
    cdef double x, y, z, sbx, sby, sbz, mass, sbmass
    cdef double distance_squared, distance, force
    cdef double dir_vec_x, dir_vec_y, dir_vec_z
    cdef double dir_vec_x_sqrd, dir_vec_y_sqrd, dir_vec_z_sqrd
    cdef double unit_vec_x, unit_vec_y, unit_vec_z

    composite_x = 0
    composite_y = 0
    composite_z = 0

    x = self.position.x
    y = self.position.y
    z = self.position.z
    mass = self.mass
    selfid = self.id

    grav_const_times_mass = grav_const * mass

    for sb_tup in spacebodies_data:
        if sb_tup[0] != selfid:

            sbmass = sb_tup[1]
            sbx = sb_tup[2]
            sby = sb_tup[3]
            sbz = sb_tup[4]

            dir_vec_x, dir_vec_y, dir_vec_z = sbx - x, sby - y, sbz - z
            dir_vec_x_sqrd = dir_vec_x ** 2
            dir_vec_y_sqrd = dir_vec_y ** 2
            dir_vec_z_sqrd = dir_vec_z ** 2

            distance_squared = dir_vec_x_sqrd + dir_vec_y_sqrd + dir_vec_z_sqrd
            distance = distance_squared ** 0.5

            force = (grav_const_times_mass * sbmass) / distance_squared

            unit_vec_x, unit_vec_y, unit_vec_z = dir_vec_x / distance, dir_vec_y / distance, dir_vec_z / distance

            composite_x += force * unit_vec_x
            composite_y += force * unit_vec_y
            composite_z += force * unit_vec_z

    return (composite_x, composite_y, composite_z)
