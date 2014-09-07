from numpy import *
from pylab import *
import chains
import matplotlib.cm as cmx
import matplotlib.colors as colors
import os, sys
import planckStyleNG

# rootparam = ['WMAP_LCDM_fixtheta', 'omegabh2', None, 1.025]
# rootparam = ['WMAP_LCDM', 'H0', None, 1.025]

# rootparam = ['LCDM_base_planck_lowl_lowLike_highL', 'H0', None]

# rootparam = ['sterile', 'nnu', [3.046, 3.3]]  # [0.27, 0.35]

# rootparam = ['LCDM_base_planck_lowl_lowLike_highL', 'omegach2', None]  # [0.27, 0.35]
# rootparam = ['LCDM_base_Alens_planck_lowl_lowLike', 'Alens', None]  # [0.9, 1.6]]
# # rootparam = ['lmax2000', 'Alens', None]

# rootparam = ['omegak_planck_lowl_lowLike_highL', 'omegak', [-0.16, 0.02]]
# rootparam = ['mnu_planck_lowl_lowLike_highL', 'mnu', [0, 0.74]]
# rootparam = ['mnu_planck_lowl_lowLike_highL', 'mnu', r'$\sum m_\nu\ [\rm{eV}]$', None]

rootparam = ['LCDM_base_planck_lowl_lowLike_highl', 'ns', None]  # [0.27, 0.35]
# rootparam = ['LCDM_base_planck_lowl_lowLike_highL', 'tau', None]  # [0.27, 0.35]

# rootparam = ['LCDM_base_planck_lowl', 'omegamh2', r'$\Omega_{\rm m} h^2$', [0.137, 0.146]]  # [0.27, 0.35]
# rootparam = ['LCDM_base_planck_lowl', 'H0', r'$H_0$', [65, 70]]  # [0.27, 0.35]

figsize = 3.5
Lmax = 2500

difference = False
export = True
colour_density = False

params = chains.loadChains(r'C:\tmp\Planck\final_nominal\cl_chains' + os.sep + rootparam[0], ignore_frac=0, separate_chains=False, no_stat=True)

nLines = 10
alpha_add = 0.6
# Lmax=410
# nLines = 800
# alpha_add = 0.02


if alpha_add < 1:
    width = 0.6 * figsize / 3.5
else: width = 1

ls = np.array(range(2, Lmax + 1))
ix = params.paramNames.numberOfName(rootparam[1])
colorlabel = '$' + params.paramNames.parWithName(rootparam[1]).label + '$'

colorsamps = params.samples[:, ix]

cmap = cm.jet


def bf_different():
    figure(figsize=(figsize, figsize * 3 / 4))
    calib = 1  # .025  # 1.025
    if len(rootparam) > 3:calib = rootparam[3]
    planck = np.loadtxt("C:\\tmp\\Planck\\final_nominal\\planck_data_points.txt")

#    bffile = 'base_Alens_planck_217onlylmax2000_lowl_lowLike'
#    bffile = 'base_planck_217onlylmax2000_lowl_lowLike'
    bffile = 'base_planck_lowl_lowLike'
#    bffile = 'base_planck_lowl_lowLike_fixtheta'
    bf = np.loadtxt("C:\\tmp\\Planck\\final_nominal\\" + bffile + '.bestfit_cl')
#    bf[3:Lmax - 1 + 3, 1] = bf[:Lmax - 1, 1]
    bf_cl = bf[:Lmax - 1, 1] * ls
    if difference:
        for i in range(5, Lmax - 5):
            if bf_cl[i] > bf_cl[i + 2] and bf_cl[i] > bf_cl[i - 2, ]:
                axvline(i, color='k', ls=':')
            if bf_cl[i] <= bf_cl[i + 2] and bf_cl[i] <= bf_cl[i - 2]:
                axvline(i, color='b', ls=':')
#    plot(bf[:Lmax, 0], bf[:Lmax, 1])
#    return

    c = chains.loadChains(r'C:\tmp\Planck\final_nominal\cl_chains' + os.sep + rootparam[0], ignore_frac=0, separate_chains=False, no_stat=True, ext='.TT')
    bf_cl = bf[:Lmax - 1, 1]
    norm = rootparam[2]
    if norm is not None:
        cNorm = colors.Normalize(vmin=norm[0], vmax=norm[1])
    else:
        cNorm = colors.Normalize(vmin=min(colorsamps[:]), vmax=max(colorsamps[:]))

    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    scalarMap.set_array(params.samples[:, ix])
    if difference: Lpow = 1
    else: Lpow = 0
    scaling = ls ** Lpow
    samps = c.samples[:, :Lmax - 1].copy()
    for i in range(nLines):
        CL = samps[i, :] / calib
        if difference: CL -= bf_cl
        samps[i, :] = CL * scaling

    if colour_density:
        CLdensity(samps, colorsamps, scalarMap, nbins=700, deltaL=5)
    else:
        for i in range(nLines):
            CL = samps[i, :]
            colorVal = scalarMap.to_rgba(colorsamps[i])
            plot(ls, samps[i, :] , linewidth=width, color=colorVal, alpha=alpha_add)

    cb = colorbar(scalarMap, norm=cNorm)
    cb.set_label(colorlabel, rotation=-90)
    axhline(0, color='k')
    xlabel(r'$\ell$')
    if difference: ylabel(r'$\ell\Delta D_\ell/\mu{\rm K}^2$')
    else: ylabel(r'$ D_\ell/\mu{\rm K}^2$')
    Lplanck = np.array([int(i) for i in planck[:, 1]])
    pts = planck[:, 4]
    if difference: pts -= bf_cl[Lplanck - 2]
#    else:
#        plot(ls, bf_cl, '-k', linewidth=1)
    plot_data(Lplanck, Lplanck ** Lpow * pts, Lplanck ** Lpow * planck[:, 5], width=figsize / 8.)
    if export:
        if difference:
            savefig(r'z://' + rootparam[0] + '_' + rootparam[1] + '_TTdiff_to_' + bffile + (('calib_' + str(calib), '')[calib == 1]) + '.pdf', bbox_inches='tight')
        else:
            savefig(r'z://' + rootparam[0] + '_' + rootparam[1] + '_TT_best' + bffile + (('calib_' + str(calib), '')[calib == 1]) + '.pdf', bbox_inches='tight')


def lowL():

    figure(figsize=(figsize, figsize * 3 / 4))
    planck = np.loadtxt("C:\\tmp\\Planck\\final_nominal\\planck_points_v6_clean.dat")

    c = chains.loadChains(r'C:\tmp\Planck\final_nominal\cl_chains' + os.sep + rootparam[0], ignore_frac=0, separate_chains=False, no_stat=True, ext='.TT')
    norm = rootparam[2]
    if norm is not None:
        cNorm = colors.Normalize(vmin=norm[0], vmax=norm[1])
    else:
        cNorm = colors.Normalize(vmin=min(colorsamps[:]), vmax=max(colorsamps[:]))

    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    scalarMap.set_array(params.samples[:, ix])
    for i in range(nLines):
        CL = c.samples[i, :Lmax - 1]  # * 1.02
        colorVal = scalarMap.to_rgba(params.samples[i, ix])
        plot(ls, CL , linewidth=width, color=colorVal, alpha=alpha_add)
    cb = colorbar(scalarMap, norm=cNorm)
    labels = [label.get_text() for label in cb.ax.yaxis.get_ticklabels()[::2]]
    cb.ax.yaxis.set_ticks(cb.ax.yaxis.get_ticklocs()[::2])
    cb.ax.yaxis.set_ticklabels(labels)

    cb.set_label(colorlabel, rotation=-90)
    for ticklabel in cb.ax.get_yticklabels():
        ticklabel.set_rotation(-90)

    xlabel(r'$\ell$')
    ylabel(r'$\mathcal{D}_\ell\,  [\mu{\rm K}^2]$')
    xlim([1, 49.8])
    ylim([0, 2500])

    meanCL = mean(c.samples[:, :Lmax - 1], 0)
    confid = np.zeros((Lmax - 1, 2))
    for i in range(Lmax - 1):
        confid[i, :] = np.percentile(c.samples[:, i], [16, 84 ])
    plot(ls, meanCL , linewidth=0.6, color='k', ls='-', zorder=2)
    dashes = [3, 1]
    P, = plot(ls, confid[:, 0] , linewidth=0.3, color='k', ls='--', zorder=2)
    P.set_dashes(dashes)
    P, = plot(ls, confid[:, 1] , linewidth=0.3, color='k', ls='--', zorder=2)
    P.set_dashes(dashes)


    col = 'k'
    planck[2, :]
    plot(planck[:, 0], planck[:, 1], 'o' , color=col, ls='None', markersize=1)
    dat = planck
    for i in range(planck.shape[0]):
        gca().add_line(Line2D((dat[i, 0], dat[i, 0]), (dat[i, 1] - dat[i, 3], dat[i, 1] + dat[i, 2]), color=col, linewidth=0.5))

    savefig(r'z://' + rootparam[0] + '_' + rootparam[1] + '_lowL.pdf', bbox_inches='tight')

def plot_data(L, d, error, width=1):
    col = 'k'
    plot(L, d, 'o' , color=col, ls='None', markersize=width * 2)
    for i in range(L.shape[0]):
        err = error[i]
        gca().add_line(Line2D((L[i], L[i]), (d[i] - err, d[i] + err), color=col, linewidth=width))


def CLdensity(samples, param, scalarMap, nbins=300, deltaL=1):
#    H, xedges, yedges = np.histogram2d(x, y, bins=(100, 20), weights=param)
    minCl = 0
    maxCl = np.max(samples) * 1.02
    delta = (maxCl - minCl) / (nbins - 1)
    nL = (Lmax - 1) / deltaL
    sums = np.zeros((nL, nbins))
    psums = np.zeros((nL, nbins))

    for i in  range(param.shape[0]):
        indices = np.floor((samples[i, ::deltaL] - minCl) / delta).astype(int)
        for j in range(nL):
            psums[j, indices[j]] += param[i]
            sums[j, indices[j]] += 1

    window = np.arange(-2.5, 2.5, 0.5)
    window = np.exp(-window ** 2 / 2)
    window /= np.sum(window)
    for i in range(nL):
        counts = sums[i, :]
        psums[i, :] = np.convolve(psums[i, :], window, 'same')
        sums[i, :] = np.convolve(counts, window, 'same')
        positive = np.where(counts > 0)
        psums[i, positive] /= counts[positive]
        sums[i, :] /= np.max(sums[i, :])

        maxix = np.max(positive)
        minix = np.min(positive)
        psums[i, maxix:] = psums[i, maxix]
        psums[i, :minix] = psums[i, minix]

    if deltaL == 1:
        for i in range(nbins):
            sums[:, i] = np.convolve(sums[:, i] , window, 'same')

    colorVal = np.zeros((nbins, nL, 4))
    for i in range(nL):
            for j in range(nbins):
                colorVal[j, i, :] = scalarMap.to_rgba(psums[i, j], sums[i, j] ** 0.25)
#    colorVal = scalarMap.to_rgba(psums, sums)
    im = imshow(colorVal, origin='lower', aspect='auto', extent=[2, Lmax, 0, maxCl])
    im.set_interpolation('bicubic')


def phi_plot():

    figure(figsize=(figsize, figsize * 3 / 4))

    dat = np.loadtxt("C:\\tmp\\Planck\\final_nominal\\lensing_points_for_AML.dat")

    c = chains.loadChains(r'C:\tmp\Planck\final_nominal\cl_chains' + os.sep + rootparam[0], ignore_frac=0, separate_chains=False, no_stat=True, ext='.PhiPhi')
    meanCL = mean(c.samples[:, :Lmax - 1], 0)
    confid = np.zeros((Lmax - 1, 2))
    for i in range(Lmax - 1):
        confid[i, :] = np.percentile(c.samples[:, i], [16, 84 ])
#    stdCL = std(c.samples[:, :Lmax - 1], 0)

    print c.numrows, c.numrows / nLines

    indices = range(0, c.numrows, max(1, c.numrows / nLines))
    norm = rootparam[2]
    if norm is not None:
        cNorm = colors.Normalize(vmin=norm[0], vmax=norm[1])
    else:
        cNorm = colors.Normalize(vmin=min(colorsamps[indices]), vmax=max(colorsamps[indices]))

    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    scalarMap.set_array(colorsamps[indices])

    if colour_density:
        CLdensity(c.samples[:, :Lmax - 1], colorsamps, scalarMap)
    else:
        for i in indices:
            CL = c.samples[i, :Lmax - 1]
            colorVal = scalarMap.to_rgba(colorsamps[i])
            plot(ls, CL , linewidth=width, color=colorVal, alpha=alpha_add, zorder=0)

    plot(ls, meanCL , linewidth=0.6, color='k', ls='-', zorder=2)
    dashes = [3, 1]
    P, = plot(ls, confid[:, 0] , linewidth=0.3, color='k', ls='--', zorder=2)
    P.set_dashes(dashes)
    P, = plot(ls, confid[:, 1] , linewidth=0.3, color='k', ls='--', zorder=2)
    P.set_dashes(dashes)

    cb = colorbar(scalarMap, norm=cNorm)
    cb.set_label(colorlabel, rotation=-90)

    # hack for errors coming out transparent
#    plot_data(dat[:,2],dat[:,3],dat[i,4]))
    col = 'k'
    plot(dat[:, 2], dat[:, 3], 'o' , color=col, ls='None', markersize=1)
    for i in range(dat.shape[0]):
        err = dat[i, 4]
        gca().add_line(Line2D((dat[i, 2], dat[i, 2]), (dat[i, 3] - err, dat[i, 3] + err), color=col, linewidth=0.5))
        gca().add_line(Line2D((dat[i, 0], dat[i, 1]), (dat[i, 3] , dat[i, 3]), color=col, linewidth=0.5))

#    errorbar(dat[:, 2], dat[:, 3], dat[:, 4], (dat[:, 0] - dat[:, 1]) / 2, color='k', marker='o', zorder=5000, ls='None', alpha=1, markersize=1)
    xlabel(r'$L$')
    ylabel(r'$[L(L+1)]^2C_L^{\phi\phi}/2\pi$')
#    setp(getp(cb.ax, 'ymajorticklabels'), fontsize=self.settings.colorbar_axes_fontsize)

    xlim([2, Lmax])
    if rootparam[1] == 'omegak': ylim([0, 2.2e-7])
    else: ylim([0, 1.8e-7])

    savefig(r'z://' + rootparam[0] + '_' + rootparam[1] + '_Cphiphi.pdf', bbox_inches='tight')


# phi_plot()
# bf_different()
lowL()
show()

# xlim([800, Lmax])
# ylim([800, Lmax])
