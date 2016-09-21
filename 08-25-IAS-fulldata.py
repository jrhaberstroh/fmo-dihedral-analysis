%pylab
import matplotlib.patheffects as path_effects
import MDAnalysis.coordinates
import os
import pandas as pd
from numpy.linalg import eigh

def unravel_dihedrals(x_t):
    ## Unravel dihedral list from (x,y,z) to (sin,cos)
    padding = 0
    even_shape=(x_t.shape[1] % 2 == 0)
    if even_shape:
        ## Check if there are two zeros padded at the end.
        ## If there are, remove one from the bulk value.
        end_two = x_t[-1, -1, -2:]
        if all(end_two == 0.0):
            padding = -1
    else:
        ## If the array shape[1] is odd, pad the array by one, and set the 
        ##   padding value to remove the last two dihedrals
        blank = np.zeros( (x_t.shape[0], 1, x_t.shape[2]) )
        x_t = np.concatenate(x_t, blank)
        padding = -2
    N_dihedral_bulk = (x_t.shape[1]) // 2 * 3
    N_dihedral = N_dihedral_bulk + padding
    x_t = x_t.reshape( (x_t.shape[0], N_dihedral_bulk, 2) )
    return(x_t[:, :N_dihedral, :])

def build_angles(x_t):
    cos = x_t[:, 0::2]
    sin = x_t[:, 1::2]
    angle = np.arctan2(sin, cos)
    return angle

def symmuncertainty(X, Y, nbins=20, xbinlim=[-np.pi, np.pi], ybinlim=[-np.pi, np.pi]):
    xbins = np.linspace(xbinlim[0], xbinlim[1], nbins)
    ybins = np.linspace(ybinlim[0], ybinlim[1], nbins)
    myhist, _, _ = np.histogram2d(X, Y, bins=[xbins, ybins])
    myhist = myhist / np.sum(myhist, axis=(0,1))
    myhist_x = myhist.sum(axis=0)
    myhist_y = myhist.sum(axis=1)
    loghist = np.log(myhist)
    loghist_x = np.log(myhist.sum(axis=0))
    loghist_y = np.log(myhist.sum(axis=1))
    I_summand = myhist * (loghist - loghist_x[:, np.newaxis] - loghist_y[np.newaxis, :])
    I = np.sum(I_summand[np.isfinite(I_summand)])
    Hxy_summand = myhist * loghist
    Hx_summand = myhist_x * loghist_x
    Hy_summand = myhist_y * loghist_y
    Hxy = - np.nansum ( Hxy_summand[np.isfinite(Hxy_summand)] )
    Hx = - np.nansum ( Hx_summand[np.isfinite(Hx_summand)] )
    Hy = - np.nansum ( Hy_summand[np.isfinite(Hy_summand)] )
    return 2 * (Hx + Hy - Hxy) / (Hx + Hy)


##########################################################################################
## Load the TRR data of dihedrals
##########################################################################################
trr50 ="/home/jhaberstroh/data/2016-06-fmo500ns/2016-08-cdc_rerun/dihedrals_50nsX10ps.trr"
trr170="/home/jhaberstroh/data/2016-06-fmo500ns/2016-08-cdc_rerun/dihedrals_170nsX10ps.trr"
n_dihedrals = 1014
dih_t = None

## Load 50ns
TRR_temp = MDAnalysis.coordinates.reader(trr50, convert_units=False)
dih_t_temp = np.zeros([TRR_temp.n_frames, TRR_temp.n_atoms, 3])
for i in xrange(TRR_temp.n_frames):
    dih_t_temp[i, :, :] = TRR_temp[i].positions[:,:]

TRR_temp.close()
dih_t_temp = unravel_dihedrals(dih_t_temp)
dih_t_temp = dih_t_temp.reshape([-1, n_dihedrals * 2])
## Convert columns to labels
dih_t_temp = pd.DataFrame(dih_t_temp)
newcols = dih_t_temp.columns.tolist()
newcols[::2] =  ["d{}cos".format(i+1) for i in xrange(n_dihedrals)]
newcols[1::2] =  ["d{}sin".format(i+1) for i in xrange(n_dihedrals)]
dih_t_temp.columns = newcols
## Add timeseries column and dataset label column as first and second columns
dih_t_temp['time, ns'] = dih_t_temp.index * .01
dih_t_temp['dataset'] = '50nsX10ps'
cols = dih_t_temp.columns.tolist()
cols = cols[-2:] + cols[:-2]
dih_t_temp = dih_t_temp[cols]
dih_t = dih_t_temp if dih_t is None else pd.concat([dih_t, dih_t_temp])


## Load 170ns
TRR_temp = MDAnalysis.coordinates.reader(trr170, convert_units=False)
dih_t_temp = np.zeros([TRR_temp.n_frames, TRR_temp.n_atoms, 3])
for i in xrange(TRR_temp.n_frames):
    dih_t_temp[i, :, :] = TRR_temp[i].positions

TRR_temp.close()
dih_t_temp = unravel_dihedrals(dih_t_temp)
dih_t_temp = dih_t_temp.reshape([-1, n_dihedrals * 2])
## Convert columns to labels
dih_t_temp = pd.DataFrame(dih_t_temp)
newcols = dih_t_temp.columns.tolist()
newcols[::2] =  ["d{}cos".format(i+1) for i in xrange(n_dihedrals)]
newcols[1::2] =  ["d{}sin".format(i+1) for i in xrange(n_dihedrals)]
dih_t_temp.columns = newcols
## Add timeseries column and dataset label column as first and second columns
dih_t_temp['time, ns'] = dih_t_temp.index * .01
dih_t_temp['dataset'] = '170nsX10ps'
cols = dih_t_temp.columns.tolist()
cols = cols[-2:] + cols[:-2]
dih_t_temp = dih_t_temp[cols]
dih_t = dih_t_temp if dih_t is None else pd.concat([dih_t, dih_t_temp])


## Load dihedrals -> resid
dih2resid_arr = np.loadtxt('output/dihedral_list_resid.ndx').astype(int)
dih2resid = pd.DataFrame()
for i in xrange(3):
    df_temp = pd.DataFrame()
    df_temp['resid'] = dih2resid_arr[:]
    if i == 0:
        df_temp['column'] = ['d{}cos'.format(i) for i in np.arange(0,len(dih2resid_arr))]
    if i == 1:
        df_temp['column'] = ['d{}sin'.format(i) for i in np.arange(0,len(dih2resid_arr))]
    if i == 2:
        df_temp['column'] = ['a{}'.format(i) for i in np.arange(0,len(dih2resid_arr))]
    dih2resid = pd.concat([dih2resid, df_temp])



## Construct residue sequence from PDB file
res_sequence = []
with open("input/4bcl.pdb") as f:
    for l in f:
        if l[0:6] == 'SEQRES':
            res_sequence = res_sequence + (l.rstrip()[19:].split())
print(res_sequence)
print(len(res_sequence))
resid2restype = pd.DataFrame(res_sequence)
resid2restype.columns = ['restype']
resid2restype['resid'] = range(1,367)
restype = [resid2restype.restype[resid2restype.resid == rid].iat[0] for rid in dih2resid.resid]
dih2resid['restype'] = restype

## Construct an angles entries, and append to an "all_t" dataframe
sincos_block = dih_t.loc[:, dih_t.columns.str.match('d\d')].as_matrix()
angles_block = build_angles(sincos_block)
angles_block = pd.DataFrame(angles_block)
angles_labels = [ 'a{}'.format(int(entry[1:-3])) for entry in dih_t.columns[dih_t.columns.str.match('d\d')][::2] ]
angles_block.columns = angles_labels
ang_t = angles_block.copy()
ang_t['time, ns'] = np.array(dih_t.loc[:, 'time, ns'])
ang_t['dataset'] = np.array(dih_t.loc[:, 'dataset'])
all_t = pd.merge(dih_t, ang_t, how='outer', on=['time, ns', 'dataset'])

##########################################################################################
## exclude bad  and redundant dihedrals
##########################################################################################

## Compute covariance and eigenanalysis
## Remove residues 120-140 that are in the loop
## (i.e. dihedrals 308-394)
#! This only works for # dihedrals fewer than 3080
exclusion1 = r'3[1-9][0-9]'
exclusion2 = r'308|309'
exclusion3 = r'39[0-4]'
exclusion_str = exclusion1 + r'|' + exclusion2 + r'|' + exclusion3
col_drops = dih_t.columns.str.contains(exclusion_str)
dih_t = dih_t.loc[:, np.logical_not(col_drops)]

## Remove resid 8-14
## (i.e. 1-15)
exclusion1 = r'd[0-9]\D'
exclusion2 = r'd1[0-5]\D'
exclusion_str = exclusion1 + r'|' + exclusion2
col_drops = dih_t.columns.str.contains(exclusion_str)
dih_t = dih_t.loc[:, np.logical_not(col_drops)]

dih2resid.loc[dih2resid.loc[:,'resid'].isin(np.arange(  8, 14+1)), 'dropped'] = True
dih2resid.loc[dih2resid.loc[:,'resid'].isin(np.arange(120,140+1)), 'dropped'] = True




#### Small symmuncertainty example
## data_columns = angles_block.columns[160:190]
## R = np.zeros([len(data_columns), len(data_columns)])
## for i, col_i in enumerate(data_columns):
##     for j, col_j in enumerate(data_columns):
##         R[i,j] = symmuncertainty(angles_block.loc[:,col_i], angles_block.loc[:, col_j])
## plt.imshow(R, interpolation="None", cmap='viridis')
## plt.xticks(range(len(data_columns)), data_columns, rotation='vertical')


## Construct SU list on the near-diagonals
N = len(angles_block.columns)
R_all = np.zeros([N, N])
for offset in xrange(1,8):
    for site in xrange(0, N-offset):
        X = angles_block.iloc[:, site]
        Y = angles_block.iloc[:, site+offset]
        red_XY = symmuncertainty(X, Y, nbins=40)
        R_all[site, site+offset] = red_XY
        R_all[site+offset, site] = red_XY
## Create a criteria for redundancy and drop all upper-triangular redundant elements
redundant = np.where(R_all > .95)
upper_triangle = redundant[1] > redundant[0]
redundant = [redundant[0][upper_triangle], redundant[1][upper_triangle]] 
if set(redundant[0]).intersection(redundant[1]) != set():
    raise ValueError("Redundant dihedrals fail to form simple pairs; some form triplets or greater")
drop_list = angles_block.columns[redundant[1]]
angles_drop_list = drop_list
dihedral_drop_list_cos = ['d{}cos'.format(element[1:]) for element in drop_list]
dihedral_drop_list_sin = ['d{}sin'.format(element[1:]) for element in drop_list]
angles_block.drop(angles_drop_list, axis=1, inplace=True)
dih_t.drop(dihedral_drop_list_cos, axis=1, inplace=True)
dih_t.drop(dihedral_drop_list_sin, axis=1, inplace=True)

## Assign angle numbers now that dihedrals have been sparsified
dih2resid['dropped'] =   dih2resid.column.isin(angles_drop_list) + \
                         dih2resid.column.isin(dihedral_drop_list_cos) + \
                         dih2resid.column.isin(dihedral_drop_list_sin)
dih2resid.loc[:, 'resangle number'] = 0
for this_resid in unique(dih2resid.resid):
    residue_match = dih2resid.resid==this_resid
    angle_match = dih2resid.loc[:, 'column'].str.match('^a\d')
    drop_match = dih2resid.dropped == False
    all_match = residue_match * angle_match * drop_match
    angle_list = np.arange(np.sum(all_match)) + 1
    dih2resid.loc[all_match, 'resangle number'] = angle_list
    dih2resid.loc[residue_match * drop_match * 
                dih2resid.loc[:,'column'].str.match('d\d*cos'), 
                             'resangle number'] = angle_list
    dih2resid.loc[residue_match * drop_match * 
                dih2resid.loc[:,'column'].str.match('d\d*sin'), 
                             'resangle number'] = angle_list

##########################################################################################
## Build "dihedral basin parameter"
##########################################################################################
## Load ideal angle database
restype2dihedral = pd.DataFrame()
with open('input/dihedral_ideal.txt') as f:
    for l in f:
        l_arr = l.split()
        restype = l_arr[0]
        ideal_angles = l_arr[2:]
        ideal_angles = ['sp2' if this_angle == '--' else this_angle for this_angle in ideal_angles]
        for i, this_angle in enumerate(ideal_angles):
            appender = {'restype': restype,
                        'resangle number': i+1,
                        'ideal angles': this_angle}
            restype2dihedral = restype2dihedral.append(appender, ignore_index=True)
restype2dihedral.loc[:, 'resangle number'] =  \
        restype2dihedral.loc[:, 'resangle number'].astype(int)
## Merge ideal angles like a boss
dih2resid = pd.merge(dih2resid, restype2dihedral, \
        on=['restype', 'resangle number'], how='left')
## Create center separatrix for ideal angles for all non-sp2 dihedrals
sp3sp3_rows = np.logical_not(dih2resid.loc[:,'ideal angles'].isin(["sp2",NaN]))
sp3sp3_rows = dih2resid.index[sp3sp3_rows]
sp3sp3_seprows = []
for row in sp3sp3_rows:
    this_ideal = dih2resid.loc[row, 'ideal angles']
    separatrix=[]
    ideal_list = [float(x) * np.pi / 180 for x in this_ideal.split(',')]
    pairs = [(ideal_list[0], ideal_list[1]),
            (ideal_list[1], ideal_list[2]),
            (ideal_list[2], ideal_list[0])]
    for pair in pairs:
        pairmin = pair[0]
        pairmax = pair[1]
        if pairmax < pairmin:
            pairmax = pairmax + (2. * np.pi)
        this_sep = (pairmax + pairmin) / 2.
        if this_sep > np.pi:
            this_sep -= np.pi * 2.
        separatrix.append(this_sep)
    sep_str = ','.join(["{:.0f}".format(separatrix_i*180./np.pi) for separatrix_i in separatrix])
    sp3sp3_seprows.append(sep_str)
dih2resid.loc[sp3sp3_rows, 'separatrix'] = sp3sp3_seprows

## Load exceptional dihedral preferences database into lists
bad_sp3=[]
good_sp2x1=[]
good_sp2x2=[]
good_sp2xN=[]
bad_sp2=[]
file_sections={ 'bad sp3':bad_sp3, 
                'good sp2 onefold':good_sp2x1, 
                'good sp2 twofold':good_sp2x2, 
                'good sp2 nfold':good_sp2xN, 
                'sp2 trash':bad_sp2}
with open("input/dihedral_ideal_exception.txt") as f:
    this_section=None
    for l in f:
        line=l.split('#')[0].strip()
        if len(line) == 0:
            continue
        if line[0] == '[':
            section_str = line[1:-1].strip()
            this_section = file_sections[section_str]
        else:
            if this_section is None:
                raise ValueError("Encountered uncommented line before section marker")
            this_section.append(line)
## Apply the modifications from the populated lists
dih2resid.loc[:,'manual separatrix'] = False
for l in bad_sp3 + bad_sp2 + good_sp2x1:
    angle_name = "a{}".format(l.split()[0])
    dih2resid.loc[dih2resid['column'] == angle_name, ['separatrix', 'manual separatrix']] = [NaN, True]
for l in good_sp2x2  + good_sp2xN:
    angle_name = "a{}".format(l.split()[0])
    sep_values = l.split()[1]
    dih2resid.loc[dih2resid['column'] == angle_name, ['separatrix', 'manual separatrix']] = [sep_values, True]

## Define the dynamical bin state
import collections
def bin_angular_trajectory(trajectory, binedges):
    binedges_queue = collections.deque(binedges)
    def generate_permutations(values_queue):
        for i in range(len(values_queue)):
            values_queue.rotate()
            yield list(values_queue)
    all_permutations = list(generate_permutations(binedges_queue))
    ## Assert that separatrix values are sorted up to a cyclic permutation
    assert(sorted(binedges) in all_permutations)
    category_t = np.zeros(len(trajectory))
    for binindex, edgepair in enumerate(zip(binedges[:], binedges[1:]+binedges[0:1])):
        print(edgepair)
        above = (trajectory > edgepair[0]).astype(bool)
        below = (trajectory <= edgepair[1]).astype(bool)
        if edgepair[0] > edgepair[1]:
            binmatches = above + below
        else:
            binmatches = above * below
        print(np.sum(above))
        print(np.sum(below))
        print(np.sum(binmatches))
        category_t[binmatches] = binindex + 1
    return category_t
for anglr in all_t.columns[all_t.columns.str.match('^a\d')]:
    trajectory = all_t[anglr].as_matrix()
    this_separatrix = dih2resid.loc[dih2resid['column'] == anglr, 'separatrix'].iloc[0]
    print(this_separatrix)
    if this_separatrix != this_separatrix:
        continue
    binedges = [float(__itr__)*np.pi/180. for __itr__ in this_separatrix.split(',')]
    states_t = bin_angular_trajectory(trajectory, binedges)
    all_t.loc[:, 'state-{}'.format(anglr)] = states_t
    



#### Create matrix of all SU pairs
## N = len(angles_block.columns)
## R_all = np.zeros([N, N])
## for offset in xrange(0, N):
##     for site in xrange(offset, N):
##         X = angles_block.iloc[:, offset]
##         Y = angles_block.iloc[:, site]
##         red_XY = symmuncertainty(X, Y, nbins=20)
##         R_all[site, offset] = red_XY
##         R_all[offset, site] = red_XY
## myhist, _, _ = np.histogram2d(angles_block.loc[:,'a201'], angles_block.loc[:,'a202'], bins=[bins, bins])
## myhist = myhist / np.sum(myhist, axis=(0,1))
## myhist_x = myhist.sum(axis=0)
## myhist_y = myhist.sum(axis=1)
## loghist = np.log(myhist)
## loghist_x = np.log(myhist.sum(axis=0))
## loghist_y = np.log(myhist.sum(axis=1))
## I = np.nansum ( myhist * (loghist - loghist_x[:, np.newaxis] - loghist_y[np.newaxis, :]) )
## Hx = - np.nansum ( myhist_x * loghist_x )
## Hy = - np.nansum ( myhist_y * loghist_y )


################################################################################
## Load the cdc data from text
################################################################################
cdcstr50  ="/home/jhaberstroh/data/2016-06-fmo500ns/2016-08-cdc_rerun/cdc371_50nsX10ps.cdc"
cdcstr170 ="/home/jhaberstroh/data/2016-06-fmo500ns/2016-08-cdc_rerun/cdc371_170nsX10ps.cdc"

cdc50  = np.genfromtxt( cdcstr50, usecols=np.arange(2,361))
cdc170 = np.genfromtxt(cdcstr170, usecols=np.arange(2,361))
cdc50  = pd.DataFrame(cdc50)
cdc170 = pd.DataFrame(cdc170)
from rrid2resid_4BCL import rrid2resid_4BCL

column_names = ["RESID{}".format(rrid2resid_4BCL(i)) for i in xrange(350)]
column_names = column_names + ['bcl{}'.format(i) for i in xrange(367, 374)] + ['solvent', 'ion']
cdc50.columns  = column_names
cdc170.columns = column_names
cdc50['time, ns'] = cdc50.index * .01
cdc50['dataset']  = '50nsX10ps'
cdc170['time, ns'] = cdc170.index * .01
cdc170['dataset'] = '170nsX10ps'

cdc = pd.concat([cdc50, cdc170])

badres = ['RESID8','RESID9','RESID10','RESID11','RESID12','RESID13','RESID14','RESID15','RESID16'] + ['RESID120','RESID121','RESID122','RESID123','RESID124','RESID125','RESID126','RESID127','RESID128','RESID129','RESID130','RESID131','RESID132','RESID133','RESID134','RESID135','RESID136','RESID137','RESID138','RESID139','RESID140']
cdc.drop(badres, axis=1, inplace=True)
all_t = all_t.merge(cdc, on=['time, ns', 'dataset'], how='left')


## Load sp2 and sp3-exception handling data
with open("input/dihedral_ideal_exception.txt") as f:
    for l in f:
        raise NotImplementedError('blah')


##############################################################################
## Direct analysis
##############################################################################
## Plot all of the circular distributions
def circular_histogram(data, radians=True, nbins=80, vlines=[]):
    bottom = 6
    max_height = 4
    theta_bins = np.linspace(-np.pi, np.pi, nbins+1, endpoint=True)
    ang_hist, _ = np.histogram(data, bins=theta_bins)
    radii = ang_hist / float(max(ang_hist)) * max_height
    width = (2*np.pi) / nbins
    ax = plt.subplot(111, polar=True)
    bars = ax.bar(theta_bins[:-1], radii, width=width, bottom=bottom)
    # Use custom colors and opacity
    for r, bar in zip(radii, bars):
        bar.set_facecolor(plt.cm.jet(r / 10.))
        bar.set_alpha(0.8)
    for line in vlines:
        theta = [line] * 2
        r_line = [bottom, bottom+max_height]
        ax.plot(theta, r_line, '--k', linewidth=2)
    ax.set_rticks([])
    ax.set_xticks(np.pi/180. * np.linspace(-135, 135, 4))
rcParams['figure.figsize']=7,7
for anglr in all_t.columns[all_t.columns.str.match("^a\d")]:
    anglr_index = dih2resid.index[dih2resid.loc[:,'column'].str.match('^'+anglr+'$')]
    this_resid    = dih2resid.loc[anglr_index, 'resid'].iat[0]
    this_restype  = dih2resid.loc[anglr_index, 'restype'].iat[0]
    this_resangle = dih2resid.loc[anglr_index, 'resangle number'].iat[0]
    this_ideal    = dih2resid.loc[anglr_index, 'ideal angles'].iat[0]
    this_separatrix = dih2resid.loc[anglr_index, 'separatrix'].iat[0]
    separatrix = []
    ## Check for nan
    if this_separatrix == this_separatrix:
        print(this_separatrix)
        separatrix=this_separatrix.split(',')
    circular_histogram(all_t.loc[:,anglr], vlines=separatrix)
    plt.title("{}: {}{} $\chi${} = {}".format(anglr, this_resid, this_restype,
            this_resangle, this_ideal))
    plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-09/angular/{}.png".format(anglr))
    plt.clf()

## Plot timeseries with color splash according to bin
for anglr in all_t.columns[all_t.columns.str.match("^a\d")]:
    anglr_index = dih2resid.index[dih2resid.loc[:,'column'].str.match('^'+anglr+'$')]
    this_resid    = dih2resid.loc[anglr_index, 'resid'].iat[0]
    this_restype  = dih2resid.loc[anglr_index, 'restype'].iat[0]
    this_resangle = dih2resid.loc[anglr_index, 'resangle number'].iat[0]
    this_ideal    = dih2resid.loc[anglr_index, 'ideal angles'].iat[0]
    anglr_index = dih2resid.index[dih2resid.loc[:,'column'].str.match('^'+anglr+'$')]
    this_separatrix = dih2resid.loc[anglr_index, 'separatrix'].iat[0]
    separatrix = []
    ## Check for nan
    if this_separatrix == this_separatrix:
        separatrix = this_separatrix.split(',')
        separatrix = [float(separatrix_i) * np.pi / 180. for separatrix_i in separatrix]
    timeseries = all_t.loc[:, anglr]
    dihedral_bin = np.zeros(len(all_t))
    if len(separatrix) > 0:
        for pair_count in xrange(3):
            in_my_range = np.zeros(len(all_t))
            s1 = separatrix[ pair_count      ]
            s2 = separatrix[(pair_count+1) %3]
            # Deal with wrap-around
            if s1 > s2:
                in_my_range = (timeseries >= s1) + (timeseries < s2)
            else:
                in_my_range = (timeseries >= s1) * (timeseries < s2)
            in_my_range = in_my_range.as_matrix()
            dihedral_bin[in_my_range] = float(pair_count) / 2.
    # marker_colors = [plt.cm.viridis(dihedral_bin_i) for dihedral_bin_i in dihedral_bin]
    plt.scatter(all_t.loc[:, 'time, ns'], timeseries * 180. / np.pi, c=dihedral_bin, cmap='viridis', linewidth=0, alpha=.8)
    plt.title("{}: {}{} $\chi${} = {}".format(anglr, this_resid, this_restype,
            this_resangle, this_ideal))
    for s in separatrix:
        plt.hlines(s * 180. / np.pi, 0, max(all_t.loc[:, 'time, ns']), linewidth=2, linestyle='--', alpha=.4)
    plt.ylim([-180, 180])
    plt.xlim([0, max(all_t.loc[:, 'time, ns'])])
    plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-09/angular/TS-{}.png".format(anglr))
    plt.clf()


##############################################################################
## Dihedral PCA
##############################################################################
## 170ns PCA + projections
dihedrals_other = dih_t.loc[dih_t.dataset == '50nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca_t = dih_t.loc[dih_t.dataset == '170nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpcaCTR_t = dpca_t - np.mean(dpca_t, axis=0)
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
v = v * np.sign(np.sum(v, axis=0))
IPR = np.sum(np.square(np.square(v)), axis = 0) * v.shape[0]
dpcaEIG_t = np.dot(dpcaCTR_t, v)
dihedrals_ctr = dihedrals_other - np.mean(dpca_t, axis=0)
dihedrals_pca = np.dot(dihedrals_ctr, v)
## Plot correlation matrix
diag = np.sqrt(np.diagonal(dpcaCovar))
dpcaCorr = dpcaCovar / diag[:, np.newaxis] / diag[np.newaxis, :]
plt.imshow(dpcaCorr, interpolation='None', cmap='BrBG')
plt.xticks(range(len(diag)), dih_t.columns[dih_t.columns.str.match('d\d')], rotation='vertical')
plt.yticks(range(len(diag)), dih_t.columns[dih_t.columns.str.match('d\d')])
## Plot a simple 2d scatter of the first two PCs
my_cmap = plt.cm.get_cmap("viridis")
c = my_cmap(np.linspace(0,1,dpcaEIG_t.shape[0]))
scatter = plt.scatter(dpcaEIG_t[:,0], dpcaEIG_t[:,1], color=c, alpha=.3, linewidth=0)
x_pos = np.mean(dpcaEIG_t[:500,0])
y_pos = np.mean(dpcaEIG_t[:500,1])
text = plt.text(x_pos, y_pos, "170nsX10ps", fontsize=15, color='white', horizontalalignment='center')
text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])
## Plot a simple 2d scatter of the first two PCs on the cross dataset
c = my_cmap(np.linspace(0,1,dihedrals_pca.shape[0]))
scatter_other = plt.scatter(dihedrals_pca[:,0], dihedrals_pca[:,1], color=c, alpha=.3, linewidth=1, edgecolor='k')
x_mean = np.mean(dihedrals_pca[:,0])
y_mean = np.mean(dihedrals_pca[:,1])
text = plt.text(x_mean, y_mean, "50nsX10ps", fontsize=15, color='white', horizontalalignment='center')
text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/scatter_pca170ns.png")


## 50ns PCA + projections
dihedrals_other = dih_t.loc[dih_t.dataset == '170nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca_t = dih_t.loc[dih_t.dataset == '50nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpcaCTR_t = dpca_t - np.mean(dpca_t, axis=0)
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
v = v * np.sign(np.sum(v, axis=0))
IPR = np.sum(np.square(np.square(v)), axis = 0) * v.shape[0]
dpcaEIG_t = np.dot(dpcaCTR_t, v)
dihedrals_ctr = dihedrals_other - np.mean(dpca_t, axis=0)
dihedrals_pca = np.dot(dihedrals_ctr, v)
## Plot a simple 2d scatter of the first two PCs
my_cmap = plt.cm.get_cmap("viridis")
c = my_cmap(np.linspace(0,1,dpcaEIG_t.shape[0]))
scatter = plt.scatter(dpcaEIG_t[:,0], dpcaEIG_t[:,1], color=c, alpha=.3, linewidth=0)
x_pos = np.mean(dpcaEIG_t[:500,0])
y_pos = np.mean(dpcaEIG_t[:500,1])
text = plt.text(x_pos, y_pos, "50nsX10ps", fontsize=15, color='white', horizontalalignment='center')
text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])
## Plot a simple 2d scatter of the first two PCs on the cross dataset
c = my_cmap(np.linspace(0,1,dihedrals_pca.shape[0]))
scatter_other = plt.scatter(dihedrals_pca[:,0], dihedrals_pca[:,1], color=c, alpha=.3, linewidth=1, edgecolor='k')
x_mean = np.mean(dihedrals_pca[:,0])
y_mean = np.mean(dihedrals_pca[:,1])
text = plt.text(x_mean, y_mean, "170nsX10ps", fontsize=15, color='white', horizontalalignment='center')
text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/scatter_pca50ns.png")



## Pooled PCA + projections + 3d scatterplot
dpca50  = dih_t.loc[dih_t.dataset == '50nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca170 = dih_t.loc[dih_t.dataset =='170nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca50  = dpca50  - np.mean(dpca50,  axis=0)
dpca170 = dpca170 - np.mean(dpca170, axis=0)
dpcaCTR_t = np.vstack([dpca50, dpca170])
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
v = v * np.sign(np.sum(v, axis=0))
IPR = np.sum(np.square(np.square(v)), axis = 0) * v.shape[0]
dpcaEIG_t = np.dot(dpcaCTR_t, v)
dpca50_EIG  = np.dot(dpca50,  v)
dpca170_EIG = np.dot(dpca170, v)
## Plot a simple 2d scatter of the first two PCs on the cross dataset
c = my_cmap(np.linspace(0,1,dpca170_EIG.shape[0]))
scatter_other = plt.scatter(dpca170_EIG[:,0], dpca170_EIG[:,1], color=c, alpha=.3)
x_pos = np.mean(dpca170_EIG[:500,0])
y_pos = np.mean(dpca170_EIG[:500,1])
text = plt.text(x_pos, y_pos, "170nsX10ps", fontsize=15, color='white', horizontalalignment='center')
text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])
## Plot a simple 2d scatter of the first two PCs
my_cmap = plt.cm.get_cmap("viridis")
c = my_cmap(np.linspace(0,1,dpca50_EIG.shape[0]))
scatter = plt.scatter(dpca50_EIG[:,0], dpca50_EIG[:,1], color=c, alpha=.3, linewidth=1, edgecolor='k')
x_mean = np.mean(dpca50_EIG[:,0])
y_mean = np.mean(dpca50_EIG[:,1])
text = plt.text(x_mean, y_mean, "50nsX10ps", fontsize=15, color='white', horizontalalignment='center')
text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/scatter_pca-pooled.png")

## 3dplot
## setup the variables
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import time
from rotanimate import rotanimate
import numpy.random as RNG
my_cmap = plt.cm.get_cmap("viridis")
c50 = my_cmap(np.linspace(0,1,dpca50_EIG.shape[0]))
c170 = my_cmap(np.linspace(0,1,dpca170_EIG.shape[0]))
lw50 = np.array([1.0] * dpca50_EIG.shape[0] )
lw170= np.array([0.0] * dpca170_EIG.shape[0] )
lw  = np.concatenate( [ lw50[::10],         lw170[::10]] )
c   = np.vstack( [       c50[::10],          c170[::10]] )
v1  = np.concatenate( [dpca50_EIG[::10,0], dpca170_EIG[::10,0]] )
v2  = np.concatenate( [dpca50_EIG[::10,1], dpca170_EIG[::10,1]] )
v3  = np.concatenate( [dpca50_EIG[::10,2], dpca170_EIG[::10,2]] )
## Create the 3d figures
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
permute = RNG.permutation(len(v1))
scatter = ax.scatter(v1[permute], v2[permute], v3[permute], c=c[permute], linewidth=lw[permute], depthshade=False, alpha=.8)
angles = np.linspace(0,360,81)[:-1]
rotanimate(ax, angles,'/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/scatter3d_pca-pooled.gif',delay=10, width=15, height=12) 




##############################################################################
## [dihedral] X [cdc values], mutual information thingy
##############################################################################

## Symmetric uncertainty by-site
all_angles  = all_t.columns[all_t.columns.str.contains('^a\d')]
all_residues= all_t.columns[all_t.columns.str.contains('RESID') | all_t.columns.str.contains('bcl')]
MU = np.zeros((len(all_angles), len(all_residues)))
for i, dih_angle in enumerate(all_angles):
    print("{} of {} rows".format(i+1, len(all_angles)))
    for j, cdc_site in enumerate(all_residues):
        X = all_t.loc[:, dih_angle].as_matrix()
        Y = all_t.loc[:, cdc_site].as_matrix()
        ## Need to catch bcl resid where timeseries is all zeros
        if min(Y) == max(Y):
            MU[i,j] = 0
        else:
            SU = symmuncertainty(X, Y, nbins=20, ybinlim=[min(Y), max(Y)])
            MU[i,j] = SU
goodrows = np.any(MU > .07, axis = 1)
goodcols = np.any(MU > .07, axis = 0)
goodrowlabel = all_angles[goodrows]
goodrowlabel = [ANG + "(RESID{})".format(dih2resid.loc[dih2resid.column==ANG, 'resid'].iloc[0]) for ANG in goodrowlabel]
goodcollabel = all_residues[goodcols]
rcParams['figure.figsize']=15,15
plt.imshow(MU[goodrows,:][:, goodcols], interpolation='None', cmap='viridis')
plt.colorbar()
plt.yticks(range(len(goodrowlabel)), goodrowlabel)
plt.xticks(range(len(goodcollabel)), goodcollabel, rotation=70)
plt.tick_params(axis='both', which='major', labelsize=9)
plt.title("SU of sparsified sc-dihedrals with site-wise cdc371, any(SU) > .07")
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/su371_majorsites.png")

## Symmetric uncertainty aggregate
cdc_total = all_t.columns[all_t.columns.str.contains('RESID') | all_t.columns.str.contains('bcl')]
cdc_total = np.hstack([all_t.loc[:,cdc_total].as_matrix(), 
                        all_t.loc[:,'solvent'].as_matrix()[:,np.newaxis], 
                        all_t.loc[:,'ion'].as_matrix()[:,np.newaxis]])
cdc_total = np.sum(cdc_total, axis=1)
cdc_total_info = np.zeros(len(all_angles))
for i, dih_angle in enumerate(all_angles):
    X = all_t.loc[:, dih_angle].as_matrix()
    Y = cdc_total
    cdc_total_info[i] = symmuncertainty(X, Y, nbins=20, ybinlim=[min(Y), max(Y)])A

rcParams['figure.figsize'] = 8, 6
plt.hist(cdc_total_info, bins=20)
plt.xlabel("Uncertainty coefficient")
plt.ylabel("Counts, dihedral angles")
plt.title("Dihedral angles binned by SU371")
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/su371_hist.png")

## Compute highest variance sidechains
cdc_df = all_t.loc[:, all_t.columns[all_t.columns.str.contains('RESID') | all_t.columns.str.contains('bcl')]]
cdc_std = np.std(cdc_df.as_matrix(), axis=0)
plt.hist(cdc_std, bins=np.linspace(0,40,41))
plt.title("Distribution of site-wise CDC $\sigma$, pooled")
plt.xlabel("$\sigma$, cm$^{-1}$")
plt.ylabel("Counts, ($N(\sigma < 1)$ = {})".format(np.sum(cdc_std < 1)))
cdc_std_N2 = np.sum( (cdc_std > 1) * (cdc_std < 2) )
plt.ylim([0,cdc_std_N2*1.02])
plt.xlim([1,40])
plt.savefig("cdc_variance.png")

## Perform change-point analysis on most significant sidechains


################################################################################
## Analyze principal components by (resid within PC)
################################################################################
dpca_t = dih_t.loc[dih_t.dataset == '170nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpcaCTR_t = dpca_t - np.mean(dpca_t, axis=0)
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
v = v * np.sign(np.sum(v, axis=0))
dpcaEIG_t = np.dot(dpcaCTR_t, v)

dpcaEIG_0 = dpcaCTR_t[:, :] * v[np.newaxis, :, 0]
dpcaEIG_0_var = np.var(dpcaEIG_0, axis=0)

dpcaEIG_0_labels = dih_t.columns[dih_t.columns.str.match('d\d')]
sort_order = np.argsort(dpcaEIG_0_var)[::-1]
dpcaEIG_0_var    = dpcaEIG_0_var[sort_order]
dpcaEIG_0_labels = dpcaEIG_0_labels[sort_order]
dpcaEIG_0 = dpcaEIG_0[:, sort_order]
## Pi chart of dpca sitewise contribution to PC1
explode = np.zeros(len(dpcaEIG_0_labels))
explode[dpcaEIG_0_var <= 1E-5] += .1
explode[dpcaEIG_0_var <= 1E-4] += .1
explode[dpcaEIG_0_var <= 1E-3] += .1
explode[dpcaEIG_0_var <= 1E-2] += .1
N1 = np.sum( (dpcaEIG_0_var > 1E-2) )
N2 = np.sum( (dpcaEIG_0_var <= 1E-2)  & (dpcaEIG_0_var > 1E-3) )
N3 = np.sum( (dpcaEIG_0_var <= 1E-3)  & (dpcaEIG_0_var > 1E-4) )
N4 = np.sum( (dpcaEIG_0_var <= 1E-4)  & (dpcaEIG_0_var > 1E-5) )
N5 = np.sum( (dpcaEIG_0_var <= 1E-5) )
N_arr = [N1, N2, N3, N4, N5]
f_arr = np.cumsum(dpcaEIG_0_var)/np.sum(dpcaEIG_0_var)
f_arr = np.array([0] + list(f_arr[:-1]))
cs=plt.cm.viridis(f_arr)
figure(figsize=(8,8))
axes(aspect=1)
wedges, texts = plt.pie(dpcaEIG_0_var * 1E10, explode=explode, wedgeprops = {'linewidth':.1}, colors = cs)
f_N = np.zeros(len(N_arr))
N_running = 0
for i, N_i in enumerate(N_arr):
    f_N[i] = (f_arr[N_running] + f_arr[N_running + N_i - 1]) * .5
    N_running += N_i
for N_i, f_i in zip(N_arr, f_N):
    txt_x = np.cos( f_i * 2. * np.pi) * .7
    txt_y = np.sin( f_i * 2. * np.pi) * .7
    text = plt.text(txt_x, txt_y, "N = {}".format(N_i), verticalalignment='center', horizontalalignment='center', color='white')
    text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])
for i in xrange(N1):
    f_i = (f_arr[i] + f_arr[i+1]) * .5
    txt_x = np.cos( f_i * 2. * np.pi) * 1.0
    txt_y = np.sin( f_i * 2. * np.pi) * 1.0
    vert_a = 'bottom' if f_i < .5 else 'top'
    hori_a = 'left' if (f_i < .25 or f_i > .75) else ('right' if (f_i < .75 and f_i > .25) else 'center')
    print(hori_a + vert_a)
    residue_i = dih2resid.loc[dih2resid.column == dpcaEIG_0_labels[i], 'resid'].iat[0]
    text = plt.text(txt_x, txt_y, str(residue_i) + dpcaEIG_0_labels[i], verticalalignment=vert_a, horizontalalignment=hori_a)
    plt.scatter(txt_x, txt_y, color='black')


## Create figures of groups of nine Scaled Eigenvector Components
myangles = []
N_approved = N1
for dih in dpcaEIG_0_labels[:N_approved]:
    angle = 'a{}'.format(dih[1:-3])
    if not angle in myangles:
        myangles.append(angle)
plt.subplot(3,1,1)
plt.plot(dpcaEIG_0[::100, :N1], linewidth=2, alpha= .6)
plt.title("Scaled eigenvector components {}-{}".format(1, N1))
plt.subplot(3,1,2)
plt.plot(dpcaEIG_0[::100, N1:(2*N1)], linewidth=2, alpha= .6)
plt.title("Scaled eigenvector components {}-{}".format(N1+1, (2*N1)))
plt.subplot(3,1,3)
plt.plot(dpcaEIG_0[::100, (2*N1):(3*N1)], linewidth=2, alpha= .6)
plt.title("Scaled eigenvector components {}-{}".format((2*N1)+1, (3*N1)))
plt.tight_layout()
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/pc170ns_PC1_top27.png")

## Create a large figure of the first nine pc scaled eigenvector components
plt.clf()
plottr = dpcaEIG_0[:, :9] - np.amax(dpcaEIG_0[:, :9], axis=0)
plottr = plottr / np.amin(plottr, axis = 0)
plt.plot( plottr[::100,:], linewidth=2, alpha= .6)
plt.title("Scaled eigenvector components {}-{}".format(1, N1))

## Plot gif sequence of relevant angles
import jh_savefig
figures = []
for i in xrange(len(myangles)):
    fig = plt.figure(i)
    residue = dih2resid.loc[dih2resid.column == myangles[i], 'resid'].iat[0]
    plt.clf()
    plt.plot(all_t.loc[::100, 'time, ns'], all_t.loc[::100, myangles[i]], 'o', alpha=.3)
    plt.ylim([-np.pi, np.pi])
    plt.title("#{}: {}, residue {}".format(i+1, myangles[i], residue))
    figures.append(fig)
jh_savefig.autoanimate("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/pc170ns_PC1_top9.gif", figures, delay=1000)


## SU_X_PC1var scatterplot
dpcaEIG_0_varangle = np.sum(dpcaEIG_0_var.reshape((-1,2)), axis=1)
plt.scatter(dpcaEIG_0_varangle, cdc_total_info)
plt.xlabel("variance, dimensionless")
plt.ylabel("symmetric uncertainty")
plt.title("relation of sin+cos scaled variance contribution in dPCA::PC1 \n vs cdc-to-angle correlation metric")
plt.gca().set_xlim(left=0)
plt.gca().set_ylim(bottom=0)
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-02/170ns-dPC1var_X_su371.png")



##
## Participation Ratios
## 
## 170ns PCA
dihedrals_other = dih_t.loc[dih_t.dataset == '50nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca_t = dih_t.loc[dih_t.dataset == '170nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpcaCTR_t = dpca_t - np.mean(dpca_t, axis=0)
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
v = v * np.sign(np.sum(v, axis=0))
IPR = np.sum(np.square(np.square(v)), axis = 0) * v.shape[0]
v_170ns_max10 = [np.argsort(v[:,i])[-10:][::-1] for i in xrange(10)]
## 50ns PCA
dihedrals_other = dih_t.loc[dih_t.dataset == '170nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca_t = dih_t.loc[dih_t.dataset == '50nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpcaCTR_t = dpca_t - np.mean(dpca_t, axis=0)
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
v = v * np.sign(np.sum(v, axis=0))
IPR = np.sum(np.square(np.square(v)), axis = 0) * v.shape[0]
v_50ns_max10 = [np.argsort(v[:,i])[-10:][::-1] for i in xrange(10)]
## Pooled PCA
dpca50  = dih_t.loc[dih_t.dataset == '50nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca170 = dih_t.loc[dih_t.dataset =='170nsX10ps', dih_t.columns.str.match('d\d')].as_matrix()
dpca50  = dpca50  - np.mean(dpca50,  axis=0)
dpca170 = dpca170 - np.mean(dpca170, axis=0)
dpcaCTR_t = np.vstack([dpca50, dpca170])
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
v = v * np.sign(np.sum(v, axis=0))
IPR = np.sum(np.square(np.square(v)), axis = 0) * v.shape[0]
v_pooled_max10 = [np.argsort(v[:,i])[-10:][::-1] for i in xrange(10)]

from functools import reduce
reduce(np.intersect1d, (v1_170ns_max10, v1_50ns_max10, v1_pooled_max10))
def dihindex2resangle(dihindex_list):
    dihindex_list = [dih_t.columns[X] for X in dihindex_list]
    dihindex_list = [dih2resid.loc[dih2resid.loc[:, 'column'] == X, ('resid', 'restype', 'resangle number', 'column')].values[0] for X in dihindex_list]
    dihindex_list = [str(X[0]) + str(X[1]) + 'X' + str(X[2]) + "({})".format(X[3][1:-3]) for X in dihindex_list]
    return dihindex_list


dih_170Xpooled = np.intersect1d(v_170ns_max10[0], v_pooled_max10[0])
dih_170Xpooled = dihindex2resangle(dih_170Xpooled)
_ = [print(X, end=",") for X in dih_170Xpooled]

dih_50Xpooled = np.intersect1d(v_50ns_max10[0], v_pooled_max10[0])
dih_50Xpooled  = dihindex2resangle(dih_50Xpooled)
_ = [print(X, end=",") for X in dih_50Xpooled]

dih_50X170 = np.intersect1d(v_50ns_max10[0], v_170ns_max10[0])
dih_50X170 = dihindex2resangle(dih_50X170)
_ = [print(X, end=",") for X in dih_50X170]


intersection_max10x10 = reduce(np.intersect1d, [reduce(np.union1d, v_pooled_max10), 
                                                reduce(np.union1d, v_50ns_max10), 
                                                reduce(np.union1d, v_170ns_max10)])
intersection_max10xN =[reduce(np.intersect1d, [reduce(np.union1d, v_pooled_max10[0:i+1]), 
                                                reduce(np.union1d, v_50ns_max10[0:i+1]), 
                                                reduce(np.union1d, v_170ns_max10[0:i+1])]) for i in xrange(10)]
intersection_max10xN = [dihindex2resangle(X) for X in intersection_max10xN]
_ = [[print(Y, end=",") for Y in X] + [print()] for X in intersection_max10xN]


################################################################################
## Large-scale inference from CPA
################################################################################

with open('input/changepoints.txt') as f:
    label_arr = []
    t0_arr    = []
    tf_arr    = []
    for l in f:
        label_arr.append(l.split(':')[0])
        times = l.split('(')[1].split(')')[0]
        t0_arr.append(float(times.split('-')[0]))
        tf_arr.append(float(times.split('>')[1]))

for l, t0, tf in zip(label_arr, t0_arr, tf_arr): 
    all_t.loc[ (all_t.dataset == '170nsX10ps' ) & 
               (all_t.loc[:,'time, ns'] > t0) &
               (all_t.loc[:,'time, ns'] < tf), 'CPgroup'] = l

def circular_histogram_pair(ang_hist_1, ang_hist_2, bins, radians=True, vlines=[]):
    nbins = len(bins) - 1
    bottom = 6
    max_height = 4
    theta_bins = bins
    ## do it once
    radii = ang_hist_1 / float(max(ang_hist_1)) * max_height
    width = (2*np.pi) / nbins
    ax = plt.subplot(111, polar=True)
    bars = ax.bar(theta_bins[:-1], radii, width=width, bottom=bottom)
    for r, bar in zip(radii, bars):
        bar.set_facecolor(plt.cm.jet(r / 10.))
        bar.set_alpha(0.9)
    ## do it again
    radii = ang_hist_2 / float(max(ang_hist_2)) * max_height
    width = (2*np.pi) / nbins
    ax = plt.subplot(111, polar=True)
    bars = ax.bar(theta_bins[:-1], radii, width=width, bottom=bottom)
    for r, bar in zip(radii, bars):
        bar.set_facecolor(plt.cm.autumn(r / 10.))
        bar.set_alpha(0.8)
    ## Plot the annotations
    for line in vlines:
        theta = [line] * 2
        r_line = [bottom, bottom+max_height]
        ax.plot(theta, r_line, '--k', linewidth=2)
    ax.set_rticks([])
    ax.set_xticks(np.pi/180. * np.linspace(-135, 135, 4))
    

rcParams['figure.figsize']=7,7
for anglr in all_t.columns[all_t.columns.str.match("^a\d")]:
    for BIN1, BIN2 in zip('ABCDEFGHI', 'BCDEFGHIJ'):
        anglr_index = dih2resid.index[dih2resid.loc[:,'column'].str.match('^'+anglr+'$')]
        this_resid    = dih2resid.loc[anglr_index, 'resid'].iat[0]
        this_restype  = dih2resid.loc[anglr_index, 'restype'].iat[0]
        this_resangle = dih2resid.loc[anglr_index, 'resangle number'].iat[0]
        this_ideal    = dih2resid.loc[anglr_index, 'ideal angles'].iat[0]
        this_separatrix = dih2resid.loc[anglr_index, 'separatrix'].iat[0]
        separatrix = []
        ## Check for nan
        if this_separatrix == this_separatrix:
            print("{}: {}{} $\chi${} = {}".format(anglr, this_resid, this_restype,
                this_resangle, this_ideal), this_separatrix)
            separatrix=np.array([float(x) for x in this_separatrix.split(',')])
            separatrix= separatrix * np.pi / 180.
        nbins = 80
        theta_bins = np.linspace(-np.pi, np.pi, nbins+1, endpoint=True)
        ang_hist_1, _ = np.histogram(all_t.loc[all_t.CPgroup == BIN1,anglr], bins=theta_bins)
        ang_hist_2, _ = np.histogram(all_t.loc[all_t.CPgroup == BIN2,anglr], bins=theta_bins)
        ang_dist_1 = ang_hist_1.astype(float)
        ang_dist_2 = ang_hist_2.astype(float)
        ang_dist_1 = ang_dist_1 / np.sum(ang_dist_1)
        ang_dist_2 = ang_dist_2 / np.sum(ang_dist_2)
        Hsq = 1 - np.sum(np.sqrt(ang_dist_1 * ang_dist_2))
        H  = np.sqrt(Hsq)
        print("H={:.3f}".format(H))
        if H > .8 or this_resid in [298, 305, 354]:
            print("\tpassed!")
            circular_histogram_pair(ang_hist_1, ang_hist_2, theta_bins)
            plt.title("{}->{}, {}{} $\chi${} = {}, H={:.3f}".format(BIN1, BIN2, this_resid, this_restype,
                    this_resangle, this_ideal, H))
            plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-16/angular_pair/{}-{}.png".format(anglr, BIN1+BIN2 ))
            plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-16/angular_pair/timeordered-{}-{}.png".format(BIN1+BIN2, anglr))
            plt.clf()



################################################################################
## Circular-Linear correlation of dihedrals with PCA
################################################################################

all_t_cos = all_t.loc[:,all_t.columns.str.match("d\d+cos")].as_matrix()
all_t_sin = all_t.loc[:,all_t.columns.str.match("d\d+sin")].as_matrix()

r_xc = np.array([np.corrcoef(cos_i, cdc_total)[1,0] for cos_i in all_t_cos.T])
r_xs = np.array([np.corrcoef(sin_i, cdc_total)[1,0] for sin_i in all_t_sin.T])
r_cs = np.array([np.corrcoef(sin_i, cos_i)[1,0] for sin_i,cos_i in zip(all_t_sin.T, all_t_cos.T)])

rsq = ( np.square(r_xc) + np.square(r_xs) - 2. * r_xc * r_xs * r_cs ) / ( 1 - np.square(r_cs) )

myhist, binedge = np.histogram(rsq, bins=40, density=True)
plt.hist(rsq, bins=binedge, normed=True)
## TODO: Label them in a pretty stack

import scipy.stats
d1 = 2.
d2 = 220.
A_alpha01 = scipy.stats.chi2.ppf((1 - .01/456.), 2.) / d2 
r_alpha01 = A_alpha01 / (1 + A_alpha01)
AX = np.linspace(np.log(scipy.stats.chi2.ppf(.01, 2.)), np.log(scipy.stats.chi2.ppf((1 - .01/456.), 2.)), 300)
AX = np.exp(AX)
## Approximation to F-distribution for d2 >> 1 s.t. d2**d2 ~ inf and d2 >> d1 s.t. d2**d1 is computable
A0 = AX / d2
rnull = A0 / (1 + A0)
# rsq =  ( A0_num / A0_den ) / ( 1 + A0_num / A0_den )
plt.plot(rnull, scipy.stats.chi2.pdf(AX, 2.) * d2 * ( 1 / (1 - rnull) + (rnull / np.square(1 - rnull)) ) )
plt.vlines(r_alpha01, 0, max(myhist), linewidth=2, linestyle='--')
plt.text( r_alpha01, max(myhist),'$p < .01$, Bonferroni corrected', verticalalignment='top', horizontalalignment='left')
plt.title("Circular-linear correlation distribution vs. null hypothesis with 1/100 sample thinning")
plt.ylabel("probability density")
plt.xlabel("correlation score")




################################################################################
## Circular Autocorrelation
################################################################################

import scipy.signal
t170_cos = all_t.loc[all_t.dataset=='170nsX10ps',all_t.columns.str.match("d\d+cos")].as_matrix()
t170_sin = all_t.loc[all_t.dataset=='170nsX10ps',all_t.columns.str.match("d\d+sin")].as_matrix()

tau_values = []
for i, (cos_t, sin_t) in enumerate(zip(t170_cos.T, t170_sin.T)):
    T = len(cos_t)
    t_x = np.linspace(0, (T+1)*.01, T, endpoint=False)

    conv_cc = scipy.signal.fftconvolve(cos_t[:], cos_t[::-1], mode='full')[T-1:]
    conv_ss = scipy.signal.fftconvolve(sin_t[:], sin_t[::-1], mode='full')[T-1:]
    conv_cs = scipy.signal.fftconvolve(cos_t[:], sin_t[::-1], mode='full')[T-1:]
    conv_sc = scipy.signal.fftconvolve(sin_t[:], cos_t[::-1], mode='full')[T-1:]

    cumsum_cc_fwd = np.cumsum(cos_t[:] * cos_t[:])[::-1]
    cumsum_ss_fwd = np.cumsum(sin_t[:] * sin_t[:])[::-1]
    cumsum_cs_fwd = np.cumsum(cos_t[:] * sin_t[:])[::-1]
    cumsum_sc_fwd = np.cumsum(sin_t[:] * cos_t[:])[::-1]
    
    cumsum_cc_bak = np.cumsum(cos_t[::-1] * cos_t[::-1])[::-1]
    cumsum_ss_bak = np.cumsum(sin_t[::-1] * sin_t[::-1])[::-1]
    cumsum_cs_bak = np.cumsum(cos_t[::-1] * sin_t[::-1])[::-1]
    cumsum_sc_bak = np.cumsum(sin_t[::-1] * cos_t[::-1])[::-1]


    numerator   = conv_cc * conv_ss - conv_cs * conv_sc
    denominator = np.sqrt(   (cumsum_cc_fwd * cumsum_ss_fwd - cumsum_cs_fwd * cumsum_sc_fwd)
                           * (cumsum_cc_bak * cumsum_ss_bak - cumsum_cs_bak * cumsum_sc_bak)  )
    corr_func = numerator/denominator

    anglr = all_t.columns[all_t.columns.str.match('^a\d+')][i]

    anglr_index = dih2resid.index[dih2resid.loc[:,'column'].str.match('^'+anglr+'$')]
    this_resid    = dih2resid.loc[anglr_index, 'resid'].iat[0]
    this_restype  = dih2resid.loc[anglr_index, 'restype'].iat[0]
    this_resangle = dih2resid.loc[anglr_index, 'resangle number'].iat[0]
    this_ideal    = dih2resid.loc[anglr_index, 'ideal angles'].iat[0]
    this_separatrix = dih2resid.loc[anglr_index, 'separatrix'].iat[0]

    zero_or_less = np.where(corr_func <= 0)[0]
    first_zero   = zero_or_less[0]

    tau = np.sum(corr_func[0:first_zero] * .01)
    tau_values.append(tau)
    plt.plot(t_x,corr_func)
    plt.title('Angular autocorrelation for {}{} $\chi${}, $\\tau = {:.1f}$ns'.format(this_restype, this_resid, this_resangle, tau))
    plt.xlabel("time, ns")
    plt.ylim([-.1, 1])

    plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-16/dihedral_autocorr/autocorr-{}.png".format(anglr))
    plt.clf()

plt.hist(np.log10(tau_values), bins='fd')
plt.xticks(np.linspace(-2,2,9), [.01, .03, .1, .3, 1, 3, 10, 30, 100])
plt.xlabel("Tau, ns")
plt.ylabel("Counts")
plt.title("Histogram of Tau values for sidechain dynamics, 170ns dataset")
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-16/tau-histogram-semilogx.png")

bottom = 10**-.2
hist, binedge = np.histogram(tau_values, bins='fd')
plt.bar(binedge[:-1], hist - bottom, binedge[1:] - binedge[:-1], bottom=bottom)
plt.semilogy()
plt.xlabel("Tau, ns")
plt.ylabel("Counts")
plt.ylim([bottom, max(hist)*1.1])
plt.title("Histogram of Tau values for sidechain dynamics, 170ns dataset")
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-16/tau-histogram-semilogy.png")


################################################################################
## Wait time distributions
################################################################################
winners = []
global_bins = np.linspace(0, 100, 100)
for state_col in all_t.columns[all_t.columns.str.match('state-a\d')]:
    print(state_col)
    T_170 = all_t.loc[all_t['dataset'] == '170nsX10ps', 'time, ns']
    state_ts_170 = all_t.loc[all_t['dataset'] == '170nsX10ps', state_col]
    state_ts_cp = diff(state_ts_170)
    persistence_list = np.split(state_ts_cp, where(state_ts_cp != 0)[0])
    ignore_list = [False] * len(persistence_list)
    for i in xrange(len(persistence_list)):
        if len(persistence_list[i]) < 1:
            ignore_list[i] = True
    persistence_list_clean = persistence_list
    ## persistence_list_clean = []
    ## TODO: Fix this bullshit
    ## this_segment = []
    ## for i in xrange(len(persistence_list)):
    ##     if i == 0:
    ##         this_segment = persistence_list[0]
    ##         continue
    ##     elif ignore_list[i-1] == True or ignore_list[i] == True:
    ##         this_segment = this_segment + persistence_list[i]
    ##     ## When you find two non-ignored entries in a row, purge the this_segment variable
    ##     else:
    ##         persistence_list_clean.append(list(this_segment))
    ##         this_segment = persistence_list[i]
    ##     ## Also purge at the end of the list
    ##     if i == len(persistence_list) - 1:
    ##         persistence_list_clean.append(list(this_segment))
    ##         this_segment = []
    persistences = np.array([float(len(segment))*.01 for segment in persistence_list_clean])
    ## Filter out all where 150ns of persistence time comes from 30ps-2ns persistences
    good_states = np.sum(persistences[(persistences < 2.) * (persistences >.02)]) < 20.
    ## Filter out all where 150ns of persistence time comes from a single state
    global_state_time = [np.sum(state_ts_170 == (__state__+1)) for __state__ in xrange(int(max(state_ts_170)))]
    global_state_time = np.array([float(ti) * .01 for ti in global_state_time])
    good_changes = all(global_state_time < 150)
    
    winner = good_states and good_changes
    if winner:
        winners.append(state_col)
        plt.subplot(2,1,1)
        plt.hist(persistences, bins=global_bins)
        plt.title(state_col)
        plt.xlabel("time, ns")
        plt.ylabel("counts")
        plt.subplot(2,1,2)
        plt.scatter(T_170,state_ts_170)
        plt.ylabel("state")
        plt.xlabel("time, ns")
        plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-23/waitdist_winners/{}.png".format(state_col))
        plt.clf()


## Plot DNA-style sequence for dihedrals during their segments
spacer_len = 50
spacer = np.zeros((len(winners), spacer_len))
blocks = [spacer]
block_ctr = []
block_names = 'ABCDEFGHIJ'
block_size = []
for block_code in block_names:
    this_block = all_t.loc[all_t.CPgroup == block_code, winners].as_matrix().T
    print('{} has {} segments'.format(block_code, this_block.shape[1]))
    block_size.append(this_block.shape[1])
    if len(block_ctr) == 0:
        block_ctr.append(block_size[-1]/2. + spacer_len)
    else:
        block_ctr.append(block_ctr[-1] + block_size[-2]/2. + spacer_len + block_size[-1]/2.)
    blocks.append(this_block)
    blocks.append(spacer)
all_blocks = np.hstack(blocks)
import matplotlib.colors
cmap = plt.cm.Set1_r
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be black
cmaplist[0] = (.0,.0,.0,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0,4,5)
ticks  = np.linspace(.5, 3.5, 4)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
plt.imshow(all_blocks.T, aspect=1./200, interpolation='none', cmap=cmap, norm=norm)
plt.yticks(block_ctr, block_names)
plt.colorbar(ticks=ticks)
plt.title("48 important dihedral state vectors as timeseries")
plt.xticks(np.arange(0, len(winners)), [l.split('a')[-1] for l in winners], rotation=90)
plt.xlabel("Dihedral angle")

## Plot all-time DNA-style sequence for dihedrals
all_t_block = all_t.loc[all_t['dataset']=='170nsX10ps', winners].as_matrix().T
all_t_block_centers = [np.mean(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
all_t_block_starts  = [min(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
all_t_block_ends    = [max(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
rcParams['figure.figsize']=10,18
import matplotlib.colors
cmap = plt.cm.Set1_r
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be black
cmaplist[0] = (.0,.0,.0,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0,4,5)
ticks  = np.linspace(.5, 3.5, 4)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
plt.imshow(all_t_block.T, aspect=1./200, interpolation='none', cmap=cmap, norm=norm)
plt.hlines(all_t_block_starts, -.5, len(winners)-.5, linewidth=2, alpha=.4)
plt.hlines(all_t_block_ends,   -.5, len(winners)-.5, linewidth=2, alpha=.4)
for Y,labelname in zip(all_t_block_centers, [letter for letter in block_names]):
    plt.text(24, Y, labelname, verticalalignment='center', horizontalalignment='center', fontsize=16, alpha=.7)
plt.colorbar(ticks=ticks)
plt.title("48 important dihedral state vectors as timeseries")
plt.xticks(np.arange(0, len(winners)), [l.split('a')[-1] for l in winners], rotation=90, fontsize=10)
plt.xlabel("Dihedral angle")


## Plot only the SUPER WINNERS
super_winners = ['state-a73','state-a197','state-a241','state-a275','state-a467','state-a518','state-a631','state-a725','state-a769']
all_t_block = all_t.loc[all_t['dataset']=='170nsX10ps', super_winners].as_matrix().T
all_t_block_centers = [np.mean(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
all_t_block_starts  = [min(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
all_t_block_ends    = [max(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
rcParams['figure.figsize']=10,18
import matplotlib.colors
cmap = plt.cm.Set1_r
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be black
cmaplist[0] = (.0,.0,.0,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0,4,5)
ticks  = np.linspace(.5, 3.5, 4)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
plt.imshow(all_t_block.T, aspect=1./1000, interpolation='none', cmap=cmap, norm=norm)
plt.hlines(all_t_block_starts, -.5, len(super_winners)-.5, linewidth=2, alpha=.4)
plt.hlines(all_t_block_ends,   -.5, len(super_winners)-.5, linewidth=2, alpha=.4)
for Y,labelname in zip(all_t_block_centers, [letter for letter in block_names]):
    plt.text(len(super_winners)/2, Y, labelname, verticalalignment='center', horizontalalignment='center', fontsize=16, alpha=.7)
plt.colorbar(ticks=ticks)
plt.title("48 important dihedral state vectors as timeseries")
plt.xticks(np.arange(0, len(super_winners)), [l.split('a')[-1] for l in super_winners], rotation=90, fontsize=10)
plt.xlabel("Dihedral angle")



## Save chunks for testing
np.savetxt("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-23/state-a229.txt", all_t.loc[all_t['dataset'] == '170nsX10ps', 'state-a229'].as_matrix().astype(int), fmt='%d')
np.savetxt("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-23/stateblock-winners.txt", all_t_block = all_t.loc[all_t['dataset']=='170nsX10ps', winners].as_matrix().astype(int))

## Use a cleaning method to eliminate <1ns transitions
all_t_block = all_t.loc[all_t['dataset']=='170nsX10ps', winners].as_matrix()
all_t_diff = np.diff(all_t_block, axis=0)
for i in xrange(all_t_diff.shape[1]):
    changepoints = np.where(all_t_diff[:,i] != 0)[0] + 1
    this_block_segmented = np.split(all_t_block[:,i], changepoints)
    print(this_block_segmented)
    break


## Use Changeponit operation to clean up data
from __future__ import print_function
from __future__ import division
import scipy.signal
import collections
def plogp(p):
    if not isinstance(p, collections.Iterable):
        raise ValueError("p is not an iterable")
    plogp_vec = np.zeros(len(p))
    nonzero = (p!=0)
    plogp_vec[nonzero] = np.log(p[nonzero])*p[nonzero]
    return plogp_vec
def changepoint(TS, plot=False):
    Ni = [TS == i+1 for i in range(3)]
    NiF = [np.cumsum(ni) for ni in Ni]
    NiB = [np.cumsum(ni[::-1])[::-1] for ni in Ni]
    NF  =  np.arange(len(TS))+1
    NB  = (np.arange(len(TS))+1)[::-1]
    PiF = [nif/NF for nif in NiF]
    PiB = [nib/NB for nib in NiB]
    LF = np.array([plogp(nif/NF) * NF for nif in NiF])
    LB = np.array([plogp(nib/NB) * NB for nib in NiB])
    costfn = np.sum(LF, axis=0) + np.sum(LB, axis=0)
    ## Normalization to likelihood ratio:
    costfn = costfn - costfn[0]
    if plot:
        plt.plot(costfn)
    changepoint = np.argmax(costfn)
    likelihood_bonus = np.amax(costfn)
    params = [pif[changepoint] for pif in PiF]
    return changepoint, likelihood_bonus, params
def binseg_simple(TS):
    peak, value = changepoint(TS)
    if isinstance(peak, collections.Iterable):
        peak = peak[0]
    if value < truncate:
        print(truncate,'>', value)
        return [TS]
    else:
        first = binseg(TS[:peak], truncate, __depth=__depth+1)
        last  = binseg(TS[peak:], truncate, __depth=__depth+1)
        return (first + last)
def binseg_AICc(TS, N_decorr = 20, plot=False):
    depth = 1
    ## initial number of parameters and change per split
    AICK  = 3
    dAICK = 4
    ## List of integers for which elements of split_list are at depth=depth
    depth_list = [0]
    ## Contiguous list of split trajectories
    split_list = [TS]
    continue_splitting = True
    while continue_splitting:
        CPF = [changepoint(split_list[segment_num], plot=plot) for segment_num in depth_list]
        lhoodrat = [x[1] for x in CPF]
        changept = [x[0] for x in CPF]
        params   = [x[2] for x in CPF]
        ## Sort from largest to smallest
        best_splitters = np.argsort(lhoodrat)[::-1]
        accepted_global_idx=[]
        accepted_changepoint=[]
        for split_idx in best_splitters:
            AICc0 = 2*(AICK+    0) + 2*(AICK+    0)*(AICK+    0+1) / \
                    (len(TS)/N_decorr - (AICK+    0) - 1)
            AICcT = 2*(AICK+dAICK) + 2*(AICK+dAICK)*(AICK+dAICK+1) / \
                    (len(TS)/N_decorr - (AICK+dAICK) - 1)
            dAICcT = AICcT - AICc0
            print('AICc Penalty: ',dAICcT,'vs',lhoodrat[split_idx])
            if dAICcT < lhoodrat[split_idx]:
                AICK += dAICK
                accepted_global_idx.append(depth_list[split_idx])
                accepted_changepoint.append(changept[split_idx])
                print("Changepoint t={}/{} at index/depth {}/{} accepted, K={}".format(changept[split_idx],len(split_list[depth_list[split_idx]]), split_idx, depth, AICK))
        ## Check if any are accepted and apply the splits
        continue_splitting = (len(accepted_global_idx) > 0)
        if continue_splitting:
            ## Sort by global time ordering... if the array is of length to sort
            if len(accepted_global_idx) > 1:
                sort_order = np.argsort(accepted_global_idx)
                accepted_global_idx  = np.array(accepted_global_idx)[sort_order]
                accepted_changepoint = np.array(accepted_changepoint)[sort_order]
            print(type(accepted_global_idx), type(accepted_changepoint))
            new_depth_list = []
            offset = 0
            for global_idx, local_changepoint in zip(accepted_global_idx, accepted_changepoint):
                segment = split_list.pop(global_idx+offset)
                seg_first = segment[ :local_changepoint ]
                seg_last  = segment[ local_changepoint: ]
                split_list.insert(global_idx+offset, seg_last)
                split_list.insert(global_idx+offset, seg_first)
                new_depth_list.append(global_idx+offset)
                new_depth_list.append(global_idx+offset+1)
                offset += 1
            depth += 1
            depth_list = new_depth_list
    return split_list

def plotseg(segment):
    startpos = 0
    for this_seg in segment:
        plt.plot(np.arange(startpos, startpos+len(this_seg)), this_seg)
        startpos += len(this_seg)

my_binseg = binseg_AICc(a229, N_decorr=100)
params = [[np.sum(segment==i)/len(segment) for i in (1,2,3)] for segment in my_binseg]
p13 = np.array([[p[0], p[2]] for p in params])
for p in params: print("{:.2f}, {:.2f}, {:.2f}".format(p[0], p[1], p[2]))
# scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(params))




all_t_block_centers = [np.mean(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
all_t_block_starts  = [min(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
all_t_block_ends    = [max(all_t.loc[all_t.CPgroup == block_code, 'time, ns'])*100 for block_code in block_names]
rcParams['figure.figsize']=10,18
import matplotlib.colors
cmap = plt.cm.Set1_r
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be black
cmaplist[0] = (.0,.0,.0,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0,4,5)
ticks  = np.linspace(.5, 3.5, 4)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
plt.imshow(all_t_block.T, aspect=1./200, interpolation='none', cmap=cmap, norm=norm)
plt.hlines(all_t_block_starts, -.5, len(winners)-.5, linewidth=2, alpha=.4)
plt.hlines(all_t_block_ends,   -.5, len(winners)-.5, linewidth=2, alpha=.4)
for Y,labelname in zip(all_t_block_centers, [letter for letter in block_names]):
    plt.text(24, Y, labelname, verticalalignment='center', horizontalalignment='center', fontsize=16, alpha=.7)
plt.colorbar(ticks=ticks)
plt.title("48 important dihedral state vectors as timeseries")
plt.xticks(np.arange(0, len(winners)), [l.split('a')[-1] for l in winners], rotation=90, fontsize=10)
plt.xlabel("Dihedral angle")

