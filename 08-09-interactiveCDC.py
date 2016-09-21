%pylab
cdc_ti = np.loadtxt("output/cdc_c5_2016-07-29.txt", skiprows=1)
cdc_total = np.sum(cdc_ti, axis=1)
import statsmodels.nonparametric.smoothers_lowess as sm_np_low
t_ns = np.array(range(cdc_ti.shape[0])) * .1
out_10ns  = sm_np_low.lowess(cdc_total, t_ns, frac=.0588)
out_100ns = sm_np_low.lowess(cdc_total, t_ns, frac=.588)
plt.plot(t_ns, cdc_total, alpha=1, color='#fb9a99', label="raw@100ps sampling")
plt.plot(out_10ns [:,0], out_10ns [:,1], linewidth=3, color="k", label="10ns smoothing")
plt.plot(out_100ns[:,0], out_100ns[:,1], linewidth=3, color="k", alpha=.5, label="100ns smoothing")
plt.xlabel("Time, ns")
plt.ylabel("Energy gap shift, cm$^{-1}$")
plt.legend()

## Load files and colorscheme
files =["output/cdc_c1_2016-07-29.txt",
        "output/cdc_c2_2016-07-29.txt",
        "output/cdc_c3_2016-07-29.txt",
        "output/cdc_c4_2016-07-29.txt",
        "output/cdc_c5_2016-07-29.txt",
        "output/cdc_c6_2016-07-29.txt",
        "output/cdc_c7_2016-07-29.txt"]
datas = [np.loadtxt(f, skiprows=1) for f in files]
colors = ["#a6cee3",
          "#1f78b4",
          "#b2df8a",
          "#33a02c",
          "#fb9a99",
          "#e31a1c",
          "#fdbf6f"]



## Plot timeseries
for i, dat in enumerate(datas):
    for j, dat in enumerate(datas):
        if j != i:
            total = np.sum(dat, axis=1)
            t_ns = np.array(range(dat.shape[0])) * .1
            plt.plot(t_ns, total, color='k', alpha=.1)
        else:
            total = np.sum(dat, axis=1)
            t_ns = np.array(range(dat.shape[0])) * .1
            plt.plot(t_ns, total, color=colors[i])
    plt.title("Timeseries for c{}".format(i+367), fontsize=18)
    plt.xlabel("Time, ns", fontsize=18)
    plt.ylabel("Gap shift, cm$^{-1}$", fontsize=18)
    plt.ylim([-300,100])
    plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-08-12/longtime/c{}_wide.png".format(i+367))
    plt.clf()



## Plot stacked histograms
Nbin = 50
cumulative_total = np.zeros(Nbin)
bins = np.linspace(-600,200,Nbin+1)
text_disp_x = [0, 25, 0, 10, 10, -10, 35]
text_disp_y = [-20, 0, -20, 0, 0, -40, 50]
for i, dat in enumerate(datas):
    total = np.sum(dat, axis=1) * 1.7
    bctr = .5* bins[1:] + .5*bins[:-1]
    hist, bins = np.histogram(total, bins=bins)
    plt.plot(        bctr+12500, hist + cumulative_total, color=colors[i], linewidth=3, label="c{}".format(i+367))
    plt.fill_between(bctr+12500, cumulative_total, hist + cumulative_total, color=colors[i], alpha=.8)
    this_max = np.argmax(hist)
    plt.text(bctr[this_max]+12500 + text_disp_x[i], hist[this_max]*.5 + cumulative_total[this_max] + text_disp_y[i], 
            "c{}".format(i+367), horizontalalignment="center", weight='bold', color='w')
    cumulative_total = cumulative_total + hist
plt.ylabel("Absorption density, a.u.")
plt.xlabel("Absorption wavenumber, cm$^{-1}$")
plt.xticks([11950, 12200, 12450, 12700])
plt.xlim([11900, 12800])


## Plot horizontal histograms
for i, dat in enumerate(datas):
    total = np.sum(dat, axis=1)
    bins = np.linspace(-300,100,31)
    bctr = .5* bins[1:] + .5*bins[:-1]
    hist, bins = np.histogram(total, bins=bins)
    plt.plot(hist, bctr,  color=colors[i], linewidth=7, label="c{}".format(i+367))
    plt.barh(bctr, hist, linewidth=0, height=(bins[1]- bins[0]), color=colors[i], alpha=.3)
plt.ylim([-300,100])


from jh_autocorr import Ct
c371 = np.sum(datas[4], axis=1)
c371 = c371[100:] *  1.23981E-4 
c371_split = np.split(c371, 4)
ct_split = np.array([Ct(x, 100, f=.5)[1] for x in c371_split])
ct_split.shape
ct_err = np.std(ct_split, axis=0, ddof=1) / np.sqrt(4)

tbar, ctbar = Ct(c371, 100, f=.5/4)
plt.errorbar(tbar / 1000., ctbar, ct_err, linewidth=2, color=colors[4])
plt.vlines(2, -1E-5, 1.5E-5, linewidth=2, alpha=.6, linestyle='--')
plt.text(2.2, 1.0E-5, "2ns")
plt.vlines(5, -1E-5, 1.5E-5, linewidth=2, alpha=.6, linestyle='--')
plt.text(5.2, 1.0E-5, "5ns")
plt.vlines(.01, -1E-5, 1.5E-5, linewidth=2, alpha=.6, linestyle='--')
plt.text(.201, 1.0E-5, "10ps")

plt.hlines(0, -.1, 20, linewidth=2, alpha=.6)
plt.xlim([-.1, 20])
plt.ylim([-.3E-5, 1.2E-5])



