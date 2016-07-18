def rrid2resid_4BCL(rrid):
    if rrid < 0:
        raise ValueError("Bad rrid in rrid2pdb: {}".format(rrid))
    elif rrid <=50:
        return rrid + 8
    elif rrid <=158:
        return rrid + 11
    elif rrid <=198:
        return rrid + 14
    elif rrid <=348:
        return rrid + 17
    # BCL chromophore indices
    elif rrid <=355:
        rel_index = rrid - 348
        return (-rel_index)
    elif rrid ==356:
        return 0
    else:
        raise ValueError("Bad rrid in rrid2pdb: {}".format(rrid))
