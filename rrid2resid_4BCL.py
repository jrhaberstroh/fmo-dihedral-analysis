def rrid2resid_4BCL(rrid):
    if rrid < 0:
        raise ValueError("Bad rrid in rrid2pdb: {}".format(rrid))
    elif rrid <=50:
        return rrid + 8
    elif rrid <=158:
        return rrid + 11
    elif rrid <=198:
        return rrid + 14
    elif rrid <=349:
        return rrid + 17
    # BCL chromophore indices
    elif rrid <=356:
        rel_index = rrid - 349
        return (-rel_index)
    else:
        raise ValueError("Bad rrid in rrid2pdb: {}".format(rrid))

__rrid2resid_4BCL__ = [rrid2resid_4BCL(x) for x in xrange(357)]

import numpy as np
def resid2rrid_4BCL(resid):
    return np.where(__rrid2resid_4BCL__ == resid)[0]
