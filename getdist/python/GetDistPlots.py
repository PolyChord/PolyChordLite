import os, paramNames
import matplotlib
matplotlib.use('Agg')
from pylab import *

"""Plotting scripts for GetDist outputs"""


class GetDistPlotSettings:
    """Default sizes, font, styles etc settings for use by plots"""


    def __init__(self, subplot_size_inch=2, fig_width_inch=None):
        # if fig_width_inch set, forces fixed size, subplot_size_inch then just determines font sizes etc
        # otherwise width as wide as neccessary to show all subplots at specified size

        self.plot_meanlikes = False
        self.shade_meanlikes = False
        self.prob_label = None
        self.prob_y_ticks = False
        # self.prob_label = 'Probability'
        self.lineM = ['-k', '-r', '-b', '-g', '-m', '-c', '-y']
        self.plot_args = None
        self.solid_colors = ['#006FED', '#E03424', 'gray', '#009966', '#000866', '#336600', '#006633' , 'm', 'r']  # '#66CC99'
        self.default_dash_styles = {'--':(3, 2)}
        self.line_labels = True
        self.x_label_rotation = 0
        self.num_shades = 80
        self.fig_width_inch = fig_width_inch  # if you want to force specific fixed width
        self.setWithSubplotSize(subplot_size_inch)
        self.progress = False
        self.tight_layout = True
        self.no_triangle_axis_labels = True
# see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
        self.colormap = cm.Blues
        self.colormap_scatter = cm.jet
        self.colorbar_rotation = None  # e.g. -90

        self.param_names_for_labels = 'clik_latex.paramnames'
        self.tick_prune = None  # 'lower' or 'upper'
        self.tight_gap_fraction = 0.13  # space between ticks and the edge

        self.legend_loc = 'best'
        self.figure_legend_loc = 'upper center'
        self.legend_frame = True
        self.figure_legend_frame = True
        self.figure_legend_ncol = 1

        self.legend_rect_border = False
        self.legend_frac_subplot_margin = 0.2
        self.legend_frac_subplot_line = 0.1
        self.legend_fontsize = None

        self.num_contours = 2
        self.solid_contour_palefactor = 0.6
        self.alpha_filled_add = 0.85
        self.alpha_factor_contour_lines = 0.5

        self.axis_marker_color = 'gray'
        self.axis_marker_ls = '--'
        self.axis_marker_lw = 0.5

    def setWithSubplotSize(self, size_inch):
        self.subplot_size_inch = size_inch
        self.lab_fontsize = 7 + 2 * self.subplot_size_inch
        self.axes_fontsize = 4 + 2 * self.subplot_size_inch
        self.legend_fontsize = self.axes_fontsize
        self.font_size = self.lab_fontsize
        self.lw1 = self.subplot_size_inch / 3.0
        self.lw_contour = self.lw1 * 0.6
        self.lw_likes = self.subplot_size_inch / 6.0
        self.scatter_size = 3
        self.colorbar_axes_fontsize = self.axes_fontsize

    def rcSizes(self, axes_fontsize=None, lab_fontsize=None, legend_fontsize=None):
        self.font_size = rcParams['font.size']
        self.legend_fontsize = legend_fontsize or rcParams['legend.fontsize']
        self.lab_fontsize = lab_fontsize  or rcParams['axes.labelsize']
        self.axes_fontsize = axes_fontsize or rcParams['xtick.labelsize']
        self.colorbar_axes_fontsize = self.axes_fontsize - 1


defaultSettings = GetDistPlotSettings()

class Density1D():
    def bounds(self): return min(self.x), max(self.x)

class Density2D():
    def bounds(self):
        return (self.x1.min(), self.x2.min()), (self.x1.max(), self.x2.max())

class paramBounds():
    def __init__(self, fileName):
        self.upper = dict()
        self.lower = dict()
        if os.path.exists(fileName):
            textFileHandle = open(fileName)
            for line in textFileHandle:
                name, low, high = [text.strip() for text in  line.split()]
                if low != 'N': self.lower[name] = float(low)
                if high != 'N': self.upper[name] = float(high)
            textFileHandle.close()

    def getUpper(self, name):
        if name in self.upper: return self.upper[name]
        else: return None

    def getLower(self, name):
        if name in self.lower: return self.lower[name]
        else: return None


class SampleAnalysisGetDist():

    def __init__(self, plot_data):
        self.plot_data = plot_data
        self.newPlot()
        self.paths = dict()

    def newPlot(self):
        self.single_samples = dict()

    def get_density_grid(self, root, param1, param2, conts=2, likes=False):
        if likes:  res = self.load_2d(root, param1, param2, '_likes')
        else: res = self.load_2d(root, param1, param2)
        if res is None: return None
        result = Density2D()
        (result.pts, result.x1, result.x2) = res
        if conts > 0: result.contours = self.load_2d(root, param1, param2, '_cont', no_axes=True)[0:conts]
        return result

    def get_density(self, root, param, likes=False):
        result = Density1D()
        pts = self.load_1d(root, param)
        if pts is None: return None
        result.x = pts[:, 0]
        result.pts = pts[:, 1]
        if (likes): result.likes = self.load_1d(root, param, '.likes')[:, 1]
        return result

    def load_single_samples(self, root):
        if not root in self.single_samples: self.single_samples[root] = loadtxt(self.plot_data_file(root) + '_single.txt')
        return self.single_samples[root]

    def paramsForRoot(self, root, labelParams=None):
        names = paramNames.paramNames(self.plot_data_file(root) + '.paramnames')
        if labelParams is not None: names.setLabelsAndDerivedFromParamNames(labelParams)
        return names

    def boundsForRoot(self, root):
        return paramBounds(self.plot_data_file(root) + '.bounds')

    def plot_data_file(self, root):
        # find first match to roots that exist in list of plot_data paths
        if os.sep in root: return root
        path = self.paths.get(root, None)
        if path is not None: return path
        for datadir in self.plot_data:
            path = datadir + os.sep + root
            if os.path.exists(path + '.paramnames'):
                self.paths[root] = path
                return path
        return self.plot_data[0] + os.sep + root

    def plot_data_file_1D(self, root, name):
        return self.plot_data_file(root) + '_p_' + name

    def plot_data_file_2D(self, root, name1, name2):
        fname = self.plot_data_file(root) + '_2D_' + name2 + '_' + name1
        if not os.path.exists(fname):
            return self.plot_data_file(root) + '_2D_' + name1 + '_' + name2, True
        else: return fname, False

    def load_1d(self, root, param, ext='.dat'):
        fname = self.plot_data_file_1D(root, param.name) + ext
        if not hasattr(param, 'plot_data'): param.plot_data = dict()
        if not fname in param.plot_data:
            if not os.path.exists(fname): param.plot_data[fname] = None
            else: param.plot_data[fname] = loadtxt(fname)
        return param.plot_data[fname]

    def load_2d(self, root, param1, param2, ext='', no_axes=False):
        fname, transpose = self.plot_data_file_2D(root, param1.name, param2.name)
        if not os.path.exists(fname + ext): return None
        pts = loadtxt(fname + ext)
        if transpose: pts = pts.transpose()
        if no_axes: return pts
        x = loadtxt(fname + '_x')
        y = loadtxt(fname + '_y')
        if transpose: return (pts, y, x)
        else: return (pts, x, y)


class GetDistPlotter():

    def __init__(self, plot_data, settings=None):
        if settings is None: self.settings = defaultSettings
        else: self.settings = settings
        if isinstance(plot_data, basestring): self.plot_data = [plot_data]
        else: self.plot_data = plot_data
        self.sampleAnalyser = SampleAnalysisGetDist(self.plot_data)
        self.newPlot()

    def newPlot(self):
        clf()
        self.extra_artists = []
        self.contours_added = []
        self.param_name_sets = dict()
        self.param_bounds_sets = dict()
        self.sampleAnalyser.newPlot()
        self.fig = None

    def get_plot_args(self, plotno, **kwargs):
        if not self.settings.plot_args is None and len(self.settings.plot_args) > plotno:
            args = self.settings.plot_args[plotno]  #
            if args is None: args = dict()
        else: args = dict()
        args.update(kwargs)
        return args

    def get_dashes_for_ls(self, ls):
        return self.settings.default_dash_styles.get(ls, None)

    def get_default_ls(self, plotno=0):
        try:
            return self.settings.lineM[plotno]
        except IndexError:
            print 'Error adding line ' + plotno + ': Add more default line stype entries to settings.lineM'
            raise

    def get_line_styles(self, plotno, **kwargs):
        args = self.get_plot_args(plotno, **kwargs)
        if not 'ls' in args: args['ls'] = self.get_default_ls(plotno)[:-1]
        if not 'dashes' in args:
            dashes = self.get_dashes_for_ls(args['ls'])
            if dashes is not None: args['dashes'] = dashes
        if not 'color' in args:
            args['color'] = self.get_default_ls(plotno)[-1]
        if not 'lw' in args: args['lw'] = self.settings.lw1
        return args

    def get_color(self, plotno, **kwargs):
        return self.get_line_styles(plotno, **kwargs)['color']

    def get_linestyle(self, plotno, **kwargs):
        return self.get_line_styles(plotno, **kwargs)['ls']

    def get_alpha2D(self, plotno, filled, **kwargs):
        args = self.get_plot_args(plotno, **kwargs)
        if filled and plotno > 0: default = self.settings.alpha_filled_add
        else: default = 1
        return args.get('alpha', default)

    def paramNamesForRoot(self, root):
        if not root in self.param_name_sets: self.param_name_sets[root] = self.sampleAnalyser.paramsForRoot(root, labelParams=self.settings.param_names_for_labels)
        return self.param_name_sets[root]

    def paramBoundsForRoot(self, root):
        if not root in self.param_bounds_sets: self.param_bounds_sets[root] = self.sampleAnalyser.boundsForRoot(root)
        return self.param_bounds_sets[root]

    def checkBounds(self, root, name, xmin, xmax):
        d = self.paramBoundsForRoot(root)
        low = d.getLower(name)
        if low is not None: xmin = max(xmin, low)
        up = d.getUpper(name)
        if up is not None: xmax = min(xmax, up)
        return xmin, xmax

    def add_1d(self, root, param, plotno=0, **kwargs):
        param = self.check_param(root, param)
        density = self.sampleAnalyser.get_density(root, param, likes=self.settings.plot_meanlikes)
        if density is None: return None;

        kwargs = self.get_line_styles(plotno, **kwargs)
        plot(density.x, density.pts, **kwargs)
        if self.settings.plot_meanlikes:
            kwargs['lw'] = self.settings.lw_likes
            plot(density.x, density.likes, **kwargs)

        return density.bounds()

    def add_2d_contours(self, root, param1, param2, plotno=0, of=None, filled=False, color=None, ls=None, cols=None,
                        alpha=None, add_legend_proxy=True, **kwargs):
        param1, param2 = self.get_param_array(root, [param1, param2])

        density = self.sampleAnalyser.get_density_grid(root, param1, param2, conts=self.settings.num_contours, likes=False)
        if density is None: return None
        if alpha is None: alpha = self.get_alpha2D(plotno, filled, **kwargs)

        if filled:
            linestyles = ['-']
            if cols is None:
                if color is None:
                    if of is not None:color = self.settings.solid_colors[of - plotno - 1]
                    else: color = self.settings.solid_colors[plotno]
                if isinstance(color, basestring):
                    cols = [matplotlib.colors.colorConverter.to_rgb(color)]
                    for _ in range(1, len(density.contours)):
                        cols = [[c * (1 - self.settings.solid_contour_palefactor) + self.settings.solid_contour_palefactor for c in cols[0]]] + cols
                else: cols = color
            levels = sorted(np.append([density.pts.max() + 1], density.contours))
            CS = contourf(density.x1, density.x2, density.pts, levels, colors=cols, alpha=alpha, **kwargs)
            if add_legend_proxy: self.contours_added.append(Rectangle((0, 0), 1, 1, fc=CS.tcolors[1][0]))
            contour(density.x1, density.x2, density.pts, levels[:1], colors=CS.tcolors[1],
                    linewidths=self.settings.lw_contour, alpha=alpha * self.settings.alpha_factor_contour_lines, **kwargs)
        else:
            if color is None: color = self.get_color(plotno, **kwargs)
            cols = [color]
            if ls is None: ls = self.get_linestyle(plotno, **kwargs)
            # not supported custom dashes here yet

            if add_legend_proxy: self.contours_added.append(Line2D([0, 1], [0, 1], color=color, ls=ls))
            linestyles = [ls]
            kwargs = self.get_plot_args(plotno, **kwargs)
            kwargs['alpha'] = alpha
            contour(density.x1, density.x2, density.pts, density.contours, colors=cols , linestyles=linestyles, linewidths=self.settings.lw_contour, **kwargs)

        return density.bounds()

    def add_2d_shading(self, root, param1, param2):
        param1, param2 = self.get_param_array(root, [param1, param2])
        density = self.sampleAnalyser.get_density_grid(root, param1, param2, conts=0, likes=self.settings.shade_meanlikes)
        if density is None: return

#        pcolor(x1, x2, shade_dat, cmap=self.settings.colormap, vmax=shade_dat.max(), vmin=shade_dat.min())
        contourf(density.x1, density.x2, density.pts, self.settings.num_shades, cmap=self.settings.colormap)
# doing contour gets rid of annoying wehite lines
        contour(density.x1, density.x2, density.pts, self.settings.num_shades, cmap=self.settings.colormap)
#        contour(cs, hold='on')

    def plot_2d(self, roots, param_pair, shaded=True, filled=False, add_legend_proxy=True, **ax_args):
        if self.fig is None: self.make_figure()
        if isinstance(roots, basestring):roots = [roots]
        param_pair = self.get_param_array(roots[0], param_pair)
        if self.settings.progress: print 'plotting: ', [param.name for param in param_pair]
        if shaded and not filled: self.add_2d_shading(roots[0], param_pair[0], param_pair[1])
        mins = None
        for i, root in enumerate(roots):
            res = self.add_2d_contours(root, param_pair[0], param_pair[1], i, of=len(roots), filled=filled, add_legend_proxy=add_legend_proxy)
            if res is None: continue
            if mins is None: mins, maxs = res
            elif not shaded:
                mins = min(res[0], mins)
                maxs = max(res[1], maxs)
        if mins is None: return
        if not 'lims' in ax_args:
            lim1 = self.checkBounds(roots[0], param_pair[0].name , mins[0], maxs[0])
            lim2 = self.checkBounds(roots[0], param_pair[1].name , mins[1], maxs[1])
            ax_args['lims'] = [lim1[0], lim1[1], lim2[0], lim2[1]]
        self.setAxes(param_pair, **ax_args)

    def add_1d_marker(self, marker, color=None, ls=None):
        self.add_x_marker(marker, color, ls)

    def add_x_marker(self, marker, color=None, ls=None, lw=None):
        if color is None: color = self.settings.axis_marker_color
        if ls is None: ls = self.settings.axis_marker_ls
        if lw is None: lw = self.settings.axis_marker_lw
        axvline(marker, ls=ls, color=color, lw=lw)

    def add_y_marker(self, marker, color=None, ls=None, lw=None):
        if color is None: color = self.settings.axis_marker_color
        if ls is None: ls = self.settings.axis_marker_ls
        if lw is None: lw = self.settings.axis_marker_lw
        axhline(marker, ls=ls, color=color, lw=lw)

    def set_locator(self, axis, x=False, prune=None):
        if x: xmin, xmax = axis.get_view_interval()
        if (x and (abs(xmax - xmin) < 0.01 or max(abs(xmin), abs(xmax)) >= 1000)):
            axis.set_major_locator(MaxNLocator(self.settings.subplot_size_inch / 2 + 3, prune=prune))
        else:
            axis.set_major_locator(MaxNLocator(self.settings.subplot_size_inch / 2 + 4, prune=prune))


    def setAxisProperties(self, axis, x, prune):
        formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        axis.set_major_formatter(formatter)
        tick_params(axis='both', which='major', labelsize=self.settings.axes_fontsize)
        if x and self.settings.x_label_rotation != 0:setp(xticks()[1], rotation=self.settings.x_label_rotation)
        self.set_locator(axis, x, prune=prune)

    def setAxes(self, params=[], lims=None, do_xlabel=True, do_ylabel=True, no_label_no_numbers=False, pos=None, prune=None):
        if lims is not None: axis(lims)
        if prune is None: prune = self.settings.tick_prune
        ax = gca()
        self.setAxisProperties(ax.xaxis, True, prune)
        if pos is not None: ax.set_position(pos)  # # set [left, bottom, width, height] for the figure
        if do_xlabel and len(params) > 0:self.set_xlabel(params[0])
        elif no_label_no_numbers: ax.set_xticklabels([])
        if len(params) > 1:
            self.setAxisProperties(ax.yaxis, False, prune)
            if do_ylabel:self.set_ylabel(params[1])
            elif no_label_no_numbers: ax.set_yticklabels([])
        return ax

    def set_xlabel(self, param):
        xlabel(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize,
                verticalalignment='baseline', labelpad=4 + self.settings.font_size)  # test_size because need a number not e.g. 'medium'

    def set_ylabel(self, param):
        ylabel(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize)

    def plot_1d(self, roots, param, marker=None, marker_color=None, label_right=False,
                no_ylabel=False, no_ytick=False, no_zero=False, **ax_args):
        if isinstance(roots, basestring):roots = [roots]
        if self.fig is None: self.make_figure()
        param = self.check_param(roots[0], param)
        xmin = None
        for i, root in enumerate(roots):
            bounds = self.add_1d(root, param, i)
            if bounds is not None:
                if xmin is None: xmin, xmax = bounds
                else:
                    xmin = min(xmin, bounds[0])
                    xmax = max(xmax, bounds[1])
        if xmin is None: Exception('No roots have parameter: ' + param.name)
        if marker is not None: self.add_x_marker(marker, marker_color)
        if not 'lims' in ax_args:
            xmin, xmax = self.checkBounds(roots[0], param.name, xmin, xmax)
            ax_args['lims'] = [xmin, xmax, 0, 1.099]
        ax = self.setAxes([param], **ax_args)

        if self.settings.prob_label is not None and not no_ylabel:
            if label_right:
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                ylabel(self.settings.prob_label)
            else: ylabel(self.settings.prob_label)
        if no_ytick or not self.settings.prob_y_ticks: ax.set_yticks([])
        elif no_ylabel: ax.set_yticklabels([])
        elif no_zero:
            ticks = ax.get_yticks()
            if ticks[-1] > 1: ticks = ticks[:-1]
            ax.set_yticks(ticks[1:])

    def make_figure(self, nplot=1, nx=None, ny=None, xstretch=1, ystretch=1):
        self.newPlot()
        if nx is None: self.plot_col = int(round(sqrt(nplot / 1.4)))
        else: self.plot_col = nx
        if ny is None: self.plot_row = (nplot + self.plot_col - 1) / self.plot_col
        else: self.plot_row = ny
        if self.settings.fig_width_inch is not None:
            self.fig = figure(figsize=(self.settings.fig_width_inch, (self.settings.fig_width_inch * self.plot_row * ystretch) / (self.plot_col * xstretch)))
        else:
            self.fig = figure(figsize=(self.settings.subplot_size_inch * self.plot_col * xstretch, self.settings.subplot_size_inch * self.plot_row * ystretch))
        return self.plot_col, self.plot_row

    def get_param_array(self, root, in_params):
        if in_params is None or len(in_params) == 0: return self.paramNamesForRoot(root).names
        else:
            if not isinstance(in_params[0], paramNames.paramInfo): return self.paramNamesForRoot(root).parsWithNames(in_params, error=True)
        return in_params

    def check_param(self, root, param):
        if not isinstance(param, paramNames.paramInfo): return self.paramNamesForRoot(root).parWithName(param, error=True)
        return param

    def param_latex_label(self, root, param):
        p = self.check_param(root, param)
        return r'$' + p.label + r'$'

    def add_legend(self, legend_labels, legend_loc=None, line_offset=0, legend_ncol=None, colored_text=False, figure=False):
            if legend_loc is None:
                if figure: legend_loc = self.settings.figure_legend_loc
                else: legend_loc = self.settings.legend_loc
            if legend_ncol is None: legend_ncol = self.settings.figure_legend_ncol
            lines = []
            if len(self.contours_added) == 0:
                for i in enumerate(legend_labels):
                    args = self.get_line_styles(i[0] + line_offset)
                    lines.append(Line2D([0, 1], [0, 1], **args))
            else: lines = self.contours_added
            args = {'ncol':legend_ncol}
            if self.settings.legend_fontsize is not None: args['prop'] = {'size':self.settings.legend_fontsize}
            if colored_text:
                args['handlelength'] = 0
                args['handletextpad'] = 0
            if figure:
#                args['frameon'] = self.settings.figure_legend_frame
                self.legend = self.fig.legend(lines, legend_labels, legend_loc, **args)
                if not self.settings.figure_legend_frame:
                    # this works with tight_layout
                    self.legend.get_frame().set_edgecolor('none')
            else:
                args['frameon'] = self.settings.legend_frame and not colored_text
                self.legend = gca().legend(lines, legend_labels, legend_loc, **args)
            if not self.settings.legend_rect_border:
                for rect in self.legend.get_patches():
                    rect.set_edgecolor(rect.get_facecolor())
            if colored_text:
                for rect, text in zip(self.legend.get_patches(), self.legend.get_texts()):
                    rect.set_visible(False)
                    text.set_color(rect.get_facecolor())
                for line, text in zip(self.legend.get_lines(), self.legend.get_texts()):
                    line.set_visible(False)
                    text.set_color(line.get_color())

            return self.legend

    def finish_plot(self, legend_labels=[], legend_loc=None, line_offset=0, legend_ncol=None, no_gap=False, no_extra_legend_space=False, no_tight=False):
        has_legend = self.settings.line_labels and len(legend_labels) > 1
        if self.settings.tight_layout and not no_tight:
            if no_gap: tight_layout(h_pad=0, w_pad=0)
            else: tight_layout()

        if has_legend:
            if legend_ncol is None: legend_ncol = self.settings.figure_legend_ncol
            if legend_loc is None: legend_loc = self.settings.figure_legend_loc
            self.extra_artists = [self.add_legend(legend_labels, legend_loc, line_offset, legend_ncol, figure=True)]
            if self.settings.tight_layout and not no_extra_legend_space:
                frac = self.settings.legend_frac_subplot_margin + (len(legend_labels) / legend_ncol) * self.settings.legend_frac_subplot_line
                if 'upper' in legend_loc: subplots_adjust(top=1 - frac / self.plot_row)
                elif 'lower' in legend_loc: subplots_adjust(bottom=frac / self.plot_row)

    def default_legend_labels(self, legend_labels, roots):
        if legend_labels is None:
            if matplotlib.rcParams['text.usetex']:
                return [x.replace('_', '{\\textunderscore}') for x in roots]
            else: return roots
        else: return legend_labels

    def plots_1d(self, roots, params=None, legend_labels=None, legend_ncol=None, nx=None,
                 paramList=None, roots_per_param=False, share_y=None, markers=None, xlims=None):
        if roots_per_param:
            params = [self.check_param(roots[i][0], param) for i, param in enumerate(params)]
        else: params = self.get_param_array(roots[0], params)
        if paramList is not None:
            wantedParams = self.paramNameListFromFile(paramList)
            params = [param for param in params if param.name in wantedParams]
        nparam = len(params)
        if share_y is None: share_y = self.settings.prob_label is not None and nparam > 1
        plot_col, plot_row = self.make_figure(nparam, nx=nx)
        plot_roots = roots
        for i, param in enumerate(params):
            subplot(plot_row, plot_col, i + 1)
            if roots_per_param: plot_roots = roots[i]
            if markers is not None and i < len(markers): marker = markers[i]
            else: marker = None
#            self.plot_1d(plot_roots, param, no_ylabel=share_y and  i % self.plot_col > 0, marker=marker, prune=(None, 'both')[share_y])
            self.plot_1d(plot_roots, param, no_ylabel=share_y and  i % self.plot_col > 0, marker=marker)
            if xlims is not None: xlim(xlims[i][0], xlims[i][1])
            if share_y: self.spaceTicks(gca().xaxis, expand=True)

        self.finish_plot(self.default_legend_labels(legend_labels, roots), legend_ncol=legend_ncol)
        if share_y: subplots_adjust(wspace=0)
        return plot_col, plot_row

    def plots_2d(self, roots, param1=None, params2=None, param_pairs=None, nx=None, legend_labels=None, legend_ncol=None, filled=False):
        pairs = []
        if param_pairs is None:
            if param1 is not None:
                param1 = self.check_param(roots[0], param1)
                params2 = self.get_param_array(roots[0], params2)
                for param in params2:
                    if param.name != param1.name: pairs.append((param1, param))
            else: raise Exception('No parameter or parameter pairs for 2D plot')
        else:
            for pair in param_pairs:
                pairs.append((self.check_param(roots[0], pair[0]), self.check_param(roots[0], pair[1])))

        plot_col, plot_row = self.make_figure(len(pairs), nx=nx)

        for i, pair in enumerate(pairs):
            subplot(plot_row, plot_col, i + 1)
            self.plot_2d(roots, param_pair=pair, filled=filled, add_legend_proxy=i == 0)

        self.finish_plot(self.default_legend_labels(legend_labels, roots), legend_ncol=legend_ncol)

        return plot_col, plot_row

    def subplot(self, x, y, **kwargs):
        return subplot(self.plot_row, self.plot_col, x + (y - 1) * self.plot_col, **kwargs)

    def plots_2d_triplets(self, root_params_triplets, nx=None, filled=False, x_lim=None):
        plot_col, plot_row = self.make_figure(len(root_params_triplets), nx=nx)
        for i, (root, param1, param2) in enumerate(root_params_triplets):
            subplot(plot_row, plot_col, i + 1)
            self.plot_2d([root], param_pair=[param1, param2], filled=filled, add_legend_proxy=i == 0)
            if x_lim is not None:xlim(x_lim)
        self.finish_plot()
        return plot_col, plot_row

    def spaceTicks(self, axis, expand=True):
            lims = axis.get_view_interval()
            tick = [tick for tick in axis.get_ticklocs() if tick > lims[0] and tick < lims[1]]
            gap_wanted = (lims[1] - lims[0]) * self.settings.tight_gap_fraction
            if expand:
                lims = [min(tick[0] - gap_wanted, lims[0]), max(tick[-1] + gap_wanted, lims[1])]
                axis.set_view_interval(lims[0], lims[1])
            else:
                if tick[0] - lims[0] < gap_wanted: tick = tick[1:]
                if lims[1] - tick[-1] < gap_wanted:tick = tick[:-1]
            axis.set_ticks(tick)
            return tick


    def triangle_plot(self, roots, in_params=None, legend_labels=None, plot_3d_with_param=None, filled_compare=False, shaded=False):
        if isinstance(roots, basestring):roots = [roots]
        params = self.get_param_array(roots[0], in_params)
        plot_col = len(params)
        if plot_3d_with_param is not None: col_param = self.check_param(roots[0], plot_3d_with_param)
        self.make_figure(nx=plot_col, ny=plot_col)
        lims = dict()
        ticks = dict()
        for i, param in enumerate(params):
            subplot(plot_col, plot_col, i * plot_col + i + 1)
            self.plot_1d(roots, param, do_xlabel=i == plot_col - 1, no_label_no_numbers=self.settings.no_triangle_axis_labels,
                         label_right=True, no_zero=True, no_ylabel=True, no_ytick=True)
            # set no_ylabel=True for now, can't see how to not screw up spacing with right-sided y label
            if self.settings.no_triangle_axis_labels: self.spaceTicks(gca().xaxis, expand=not shaded)
            lims[i] = xlim()
            ticks[i] = gca().get_xticks()

        for i, param in enumerate(params):
            for i2 in range(i + 1, len(params)):
                param2 = params[i2]
                subplot(plot_col, plot_col, i2 * plot_col + i + 1)
                if plot_3d_with_param is not None:
                    self.plot_3d(roots, [param, param2, col_param], color_bar=False, line_offset=1,
                      do_xlabel=i2 == plot_col - 1, do_ylabel=i == 0, filled=filled_compare, no_label_no_numbers=self.settings.no_triangle_axis_labels)
                else:
                    self.plot_2d(roots, param_pair=[param, param2], do_xlabel=i2 == plot_col - 1, do_ylabel=i == 0,
                                        no_label_no_numbers=self.settings.no_triangle_axis_labels, shaded=shaded, add_legend_proxy=i == 1 and i2 == i + 1)
                gca().set_xticks(ticks[i])
                gca().set_yticks(ticks[i2])
                xlim(lims[i])
                ylim(lims[i2])

        if self.settings.no_triangle_axis_labels:subplots_adjust(wspace=0, hspace=0)
        if plot_3d_with_param is not None:
            cb = self.fig.colorbar(self.last_scatter, cax=self.fig.add_axes([0.9, 0.5, 0.03, 0.35]))
            self.add_colorbar_label(cb, col_param)

        self.finish_plot(self.default_legend_labels(legend_labels, roots),
                         legend_loc=None, no_gap=self.settings.no_triangle_axis_labels, no_extra_legend_space=True)

    def rectangle_plot(self, xparams, yparams, yroots, filled=True, ymarkers=None, xmarkers=None, **kwargs):
            self.make_figure(nx=len(xparams), ny=len(yparams))
#            f, plots = subplots(len(yparams), len(xparams), sharex='col', sharey='row')
            sharey = None
            yshares = []
            for x, xparam in enumerate(xparams):
                sharex = None
                for y, (yparam, roots) in enumerate(zip(yparams, yroots)):
#                    f.sca(plots[y, x])
                    if x > 0: sharey = yshares[y]
                    ax = self.subplot(x + 1, y + 1, sharex=sharex, sharey=sharey)
                    if y == 0: sharex = ax
                    self.plot_2d(roots, param_pair=[xparam, yparam], filled=filled, do_xlabel=y == len(yparams) - 1,
                                 do_ylabel=x == 0, add_legend_proxy=x == 0 and y == 0)
                    if ymarkers is not None and ymarkers[y] is not None: self.add_y_marker(ymarkers[y], **kwargs)
                    if xmarkers is not None and xmarkers[x] is not None: self.add_x_marker(xmarkers[x], **kwargs)
                    if y == 0: lims = xlim()
                    else: lims = (min(xlim()[0], lims[0]), max(xlim()[1], lims[1]))
                    if y != len(yparams) - 1: setp(ax.get_xticklabels(), visible=False)
                    if x != 0: setp(ax.get_yticklabels(), visible=False)
                    if x == 0: yshares.append(ax)

                sharex.set_xlim(lims)
                self.spaceTicks(sharex.xaxis)
                sharex.set_xlim(sharex.xaxis.get_view_interval())
            for ax in yshares:
                self.spaceTicks(ax.yaxis)
                ax.set_ylim(ax.yaxis.get_view_interval())
            subplots_adjust(wspace=0, hspace=0)
            self.finish_plot(no_gap=True)

    def rotate_yticklabels(self, ax=None, rotation=-90):
        if ax is None: ax = gca()
        for ticklabel in ax.get_yticklabels():
            ticklabel.set_rotation(rotation)

    def add_colorbar(self, param, orientation='vertical'):
        cb = colorbar(orientation=orientation)
        self.add_colorbar_label(cb, param)
        if self.settings.colorbar_rotation is not None:
            self.rotate_yticklabels(cb.ax, self.settings.colorbar_rotation)
            labels = [label.get_text() for label in cb.ax.yaxis.get_ticklabels()[::2]]
            cb.ax.yaxis.set_ticks(cb.ax.yaxis.get_ticklocs()[::2])
            cb.ax.yaxis.set_ticklabels(labels)
        return cb

    def add_line(self, P1, P2, zorder=0, color=None, ls=None, **kwargs):
            if color is None: color = self.settings.axis_marker_color
            if ls is None: ls = self.settings.axis_marker_ls
            gca().add_line(Line2D(P1, P2, color=color, ls=ls, zorder=zorder, **kwargs))

    def add_colorbar_label(self, cb, param):
        cb.set_label(r'$' + param.label + '$', fontsize=self.settings.lab_fontsize, rotation=-90)
        setp(getp(cb.ax, 'ymajorticklabels'), fontsize=self.settings.colorbar_axes_fontsize)

    def add_3d_scatter(self, root, in_params, color_bar=True):
        params = self.get_param_array(root, in_params)
        pts = self.sampleAnalyser.load_single_samples(root)
        names = self.paramNamesForRoot(root)
        cols = []
        for param in params:
            cols.append(names.numberOfName(param.name) + 2)
        self.last_scatter = scatter(pts[:, cols[0]], pts[:, cols[1]], edgecolors='none',
                s=self.settings.scatter_size, c=pts[:, cols[2]], cmap=self.settings.colormap_scatter)
        if color_bar: self.last_colorbar = self.add_colorbar(params[2])
        mins = [min(pts[:, cols[0]]), min(pts[:, cols[1]])]
        maxs = [max(pts[:, cols[0]]), max(pts[:, cols[1]])]
        for i in range(2):
            r = maxs[i] - mins[i]
            mins[i] -= r / 20
            maxs[i] += r / 20

        return [mins, maxs]

    def plot_3d(self, roots, in_params, color_bar=True, line_offset=0, filled=False, **ax_args):
        if isinstance(roots, basestring):roots = [roots]
        params = self.get_param_array(roots[0], in_params)
        if self.fig is None: self.make_figure()
        mins, maxs = self.add_3d_scatter(roots[0], params, color_bar=color_bar)
        for i, root in enumerate(roots[1:]):
            res = self.add_2d_contours(root, params[0], params[1], i + line_offset, filled=filled, add_legend_proxy=False)
            if res is None: continue
            mins = min(res[0], mins)
            maxs = max(res[1], maxs)
        if not 'lims' in ax_args:
            lim1 = self.checkBounds(roots[0], params[0].name , mins[0], maxs[0])
            lim2 = self.checkBounds(roots[0], params[1].name , mins[1], maxs[1])
            ax_args['lims'] = [lim1[0], lim1[1], lim2[0], lim2[1]]
        self.setAxes(params, **ax_args)

    def plots_3d(self, roots, param_sets, nx=None, filled_compare=False, legend_labels=None):
        sets = [[self.check_param(roots[0], param) for param in param_group] for param_group in param_sets]
        plot_col, plot_row = self.make_figure(len(sets), nx=nx, xstretch=1.3)

        for i, triplet in enumerate(sets):
            subplot(plot_row, plot_col, i + 1)
            self.plot_3d(roots, triplet, filled=filled_compare)
        self.finish_plot(self.default_legend_labels(legend_labels, roots[1:]), no_tight=True)
        return plot_col, plot_row

    def export(self, fname):
        savefig(fname, bbox_extra_artists=self.extra_artists, bbox_inches='tight')

    def paramNameListFromFile(self, fname):
        p = paramNames.paramNames(fname)
        return [name.name for name in p.names]



def sample_plots():
    g = GetDistPlotter('main/plot_data')
    g.settings.setWithSubplotSize(3)
    g.settings.param_names_for_labels = 'clik_latex.paramnames'

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_post_lensing', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO']
    params = g.get_param_array(roots[0], ['omegam', 'omegal', 'H0'])
    g.add_3d_scatter(roots[0], params)
    g.add_line([1, 0], [0, 1])
    g.add_2d_contours(roots[1], params[0], params[1], filled=False, color='g')
    g.add_2d_contours(roots[2], params[0], params[1], filled=True, color='#CC1100')
    g.setAxes(params, lims=[0, 1, 0, 1])
    g.export('omegam-omegal-H0_planck.pdf')

    g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'omegal', 'omegak', ], plot_3d_with_param='H0', filled_compare=True,
                    legend_labels=['Planck', 'Planck+BAO'])
    g.export(roots[0] + '.pdf')

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_post_lensing', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO', 'base_planck_CAMspec_lowl_lowLike']
    g.plots_1d(roots)
    g.export(roots[0] + '_1d.pdf')

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO']
    g.plots_2d(roots, param_pairs=[('omegabh2', 'ns'), ('logA', 'tau')], nx=2)
    g.export(roots[0] + '_2d.pdf')

    roots = ['base_omegak_planck_CAMspec_lowl_lowLike', 'base_omegak_planck_CAMspec_lowl_lowLike_post_lensing', 'base_omegak_planck_CAMspec_lowl_lowLike_BAO', 'base_planck_CAMspec_lowl_lowLike']
    g.triangle_plot(roots, ['omegabh2', 'omegach2', 'ns', 'omegal', 'omegak'], plot_3d_with_param='H0',
                    legend_labels=['Planck', 'Planck+lensing', 'Planck+BAO', 'Planck (flat)'])
    g.export(roots[0] + '_tri.pdf')

    roots = ['base_planck_CAMspec_lowl_lowLike', 'base_planck_CAMspec_lowl_lowLike_post_BAO']
    g.plots_3d(roots, [['omegabh2', 'omegach2', 'ns'], ['omegach2', 'logA', 'tau']])
    g.export(roots[0] + '_3d.pdf')


# sample_plots()
