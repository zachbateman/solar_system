'''
Pulling out expensive calculations to CYTHON
'''
import cython
# import copy

# needs to be "def" b/c the sorting lambda at the bottom doesn't work with "cpdef"!!!
# For decorator usage see:  https://stackoverflow.com/questions/4709285/cython-float-division-pyexc-zerodivisionerror-checking?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
@cython.cdivision(True)
cpdef calc_collisions(list spacebodies, list spacebodies_copy, double TOTAL_MASS, double MAX_DISTANCE_OF_IMPACT):

    cdef list spacebodies_to_combine, sb_ids_to_combine, id_list
    cdef str sbid, sb2id
    cdef double sbx, sby, sbz
    cdef double sb2x, sb2y, sb2z
    cdef double sbmass, sb2mass
    cdef double distance, mass_fraction, dist_of_impact
    cdef int i, num_sb_ids_to_combine
    cdef str id
    cdef bint added_to_sb_ids_to_combine  # "bint" is for boolean

    spacebodies_to_combine = []
    sb_ids_to_combine = []
    num_sb_ids_to_combine = 0

    for sb in spacebodies:
        sbid = sb.id
        sbx, sby, sbz = sb.position.x, sb.position.y, sb.position.z
        sbmass = sb.mass
        for sb2 in spacebodies:
            sb2id = sb2.id
            if sbid != sb2id:
                sb2x, sb2y, sb2z = sb2.position.x, sb2.position.y, sb2.position.z
                sb2mass = sb2.mass

                distance = ((sbx - sb2x) ** 2 + (sby - sb2y) ** 2 + (sbz - sb2z) ** 2) ** 0.5
                mass_fraction = (sbmass + sb2mass) / TOTAL_MASS
                dist_of_impact = (mass_fraction ** (1 / 3)) * MAX_DISTANCE_OF_IMPACT  # taking cube-root to adjust for spherical nature

                if distance < dist_of_impact:
                    if num_sb_ids_to_combine == 0:
                        sb_ids_to_combine.append([sbid, sb2id])
                        num_sb_ids_to_combine += 1
                    else:
                        added_to_sb_ids_to_combine = False
                        for i in range(num_sb_ids_to_combine):
                            if sbid in sb_ids_to_combine[i]:
                                sb_ids_to_combine[i] = sb_ids_to_combine[i] + [sb2id]
                                added_to_sb_ids_to_combine = True
                            elif sb2id in sb_ids_to_combine[i]:
                                sb_ids_to_combine[i] = sb_ids_to_combine[i] + [sbid]
                                added_to_sb_ids_to_combine = True
                        if not added_to_sb_ids_to_combine:
                            sb_ids_to_combine.append([sbid, sb2id])
                            num_sb_ids_to_combine += 1


    for id_list in sb_ids_to_combine:
        spacebodies_to_combine = [sb for sb in spacebodies_copy if sb.id in id_list]
        # DON'T USE DICTIONARY LOOKUP METHOD!!! CAUSES ISSUES!!!! (different/same object issues)
        # spacebodies_to_combine = [spacebodies_copy_dict[id] for id in id_list]
        result = spacebodies_to_combine[0]
        for sb in spacebodies_to_combine[1:]:
            result += sb
        spacebodies.append(result)

    spacebodies = [sb for sb in spacebodies if sb.id not in [id for id_list in sb_ids_to_combine for id in id_list]]
    return spacebodies
