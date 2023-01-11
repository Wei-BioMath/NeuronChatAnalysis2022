import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial','Roboto']

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)

def arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)

def arrow_path(xsupport, ysupport, t = 100, push=.8, arrows= 1):
    # Create arrow data
    arset = list(np.linspace(0, 1, arrows + 2))
    c = zip([xsupport[int(t * a * push)] for a in arset[1:-1]],
            [ysupport[int(t * a * push)] for a in arset[1:-1]])
    dt = zip([xsupport[int(t * a * push) + 1] - xsupport[int(t * a * push)] for a in arset[1:-1]],
             [ysupport[int(t * a * push) + 1] - ysupport[int(t * a * push)] for a in arset[1:-1]])
    arrowpath = zip(c, dt)
    return arrowpath

def get_plane_xyz(w, h, xs, ys, scale = 1.0, ymin_s=.1, ymax_s=.1, xmin_s = .05, xmax_s= .05):
    xs = xs * scale
    ys = ys * scale
    xdiff = max(xs) - min(xs)
    ydiff = max(ys) - min(ys)
    ymin = min(ys) - ydiff * ymin_s
    ymax = max(ys) + ydiff * ymax_s
    xmin = min(xs) - xdiff * xmin_s * (w / h)
    xmax = max(xs) + xdiff * xmax_s * (w / h)
    xx, yy = np.meshgrid([xmin, xmax], [ymin, ymax])
    zz = np.zeros(xx.shape)
    return xx, yy, zz, xdiff, ydiff

def plot_surface(ax, img, xs, ys, stride, xlim_ext=.1, ylim_ext=.0, zorder=0):
    img_w, img_h, _ = img.shape
    xdiff = max(xs) - min(xs)
    ydiff = max(ys) - min(ys)
    min_x, max_x = min(xs) - xdiff * xlim_ext, max(xs) + xdiff * xlim_ext
    min_y, max_y = min(ys) - ydiff * ylim_ext, max(ys) + ydiff * ylim_ext
    x = np.arange(min_x, max_x, (max_x - min_x) / img.shape[1])
    y = np.arange(min_y, max_y, (max_y - min_y) / img.shape[0])
    x, y = np.meshgrid(x, y)

    ax.plot_surface(x, y, np.zeros(x.shape), rstride=stride, cstride=stride, facecolors=img, zorder=zorder, antialiased=True)

def draw_curve(p1, p2, t = 100, d=.35, random_arrow_direction=False, direction=None):
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]

    direction = ((float(np.random.choice([-1, 1], 1, p =[0.5, 0.5]))) if random_arrow_direction else (1 if x1 > x2 else -1)) if not direction else direction

    midx = (x1 + x2) / 2.
    midy = (y1 + y2) / 2.

    m = (y2 - y1) / (x2 - x1)
    xx = midx + direction * math.sqrt((d*d) / (1. + 1. / (m*m)))

    yy = -(1./m) * (xx - midx) + midy
    f = np.poly1d(np.polyfit((p1[0], xx, p2[0]), (p1[1], yy, p2[1]), 2))

    xsupport = np.linspace(p1[0], p2[0], t)
    ysupport = f(xsupport)

    return xsupport, ysupport

def plot_3d_spatial_ccc(args,
                        img,
                        x2s,
                        y2s,
                        w=16,
                        h=4,
                        view_angle=-75,
                        view_height=14,
                        zoom_level= 8,
                        marker_sz=8,
                        lw = .7,
                        title_font_sz=7,
                        legend_font_sz= 6,
                        annotation_img_stride=1,
                        surface_color='#f8f9fa',
                        title_color= 'black',
                        surface_xl_ext = .05,
                        surface_xr_ext = .05,
                        surface_alpha=.18,
                        scatter_alpha= 1.,
                        box_to_anchor_x = .9,
                        box_to_anchor_y = .45,
                        box_width_scale = 1.,
                        legend_handle_sz = 40,
                        ylabels = ("Raw Image", "Segmentation", "Communication"),
                        y_label_horizontal_offset= 0,
                        y_label_vertical_offset= 0,
                        xlim_l_ext_offset= .4,
                        xlim_r_ext_offset= .4,
                        surface_z_offsets= (.0, .4, .8),
                        ylabel_z_offsets = (.02, .02, .02),
                        cell_type_in_str= True,
                        ccc_scatter_sz=100,
                        ccc_scatter_alpha=.8,
                        ccc_center_x= 100.,
                        ccc_center_y= 400.,
                        ccc_radius = 500.,
                        ccc_scatter_y_off=0,
                        ccc_d=.2,
                        ccc_t=100,
                        ccc_push= 1.8,
                        ccc_self_loop_push= 1.8,
                        mutation_scale=3,
                        zlim= (0, 1.1),
                        random_arrow_direction=False,
                        arrow_exception=False,
                        random_seed=0,
                        return_ax=False,
                        ):
    meta_fp = f"data/{args.dataset}/{args.meta_fn}"
    if not os.path.exists(meta_fp):
        print(f"meta file of the dataset not exist at: {meta_fp}, please double check the `meta_fn` settings in config.py or your input")
        return -1

    df = pd.read_csv(meta_fp)
    fig, ax = plt.subplots(1, 1,
                           figsize=(w, h),
                           subplot_kw={'projection': '3d'})
    xs = df["sdimx"].to_numpy().astype(float)
    ys = df["sdimy"].to_numpy().astype(float)

    xx, yy, zz, xdiff, ydiff = get_plane_xyz(w, h, xs, ys, xmin_s=surface_xl_ext, xmax_s=surface_xr_ext)

    zi = 0 # Layer 1
    ax.plot_surface(xx, yy, zz + surface_z_offsets[zi], color=surface_color, alpha=surface_alpha, zorder=zi)

    plot_surface(ax, img, xs,ys, annotation_img_stride, zorder=0)

    ax.text(y_label_horizontal_offset, y_label_vertical_offset, ylabel_z_offsets[zi] + surface_z_offsets[zi], ylabels[zi], 'y',
            color=title_color, fontsize=title_font_sz, zorder=1e5, ha='left', va='center')

    zi = 1 # Layer 2
    cell_types = df["cell_type"].to_numpy()
    cell_types = cell_types.astype(str) if cell_type_in_str else cell_types.astype(int)
    colormap = plt.get_cmap(args.colormap)

    uniq_cell_types = np.unique(cell_types)
    cell_type_colors = {cell_type: colormap((cid) / float(len(uniq_cell_types) - 1))
                        for cid, cell_type in
                        enumerate(uniq_cell_types)}
    zs = np.zeros_like(xs)
    for cid, cell_type in enumerate(uniq_cell_types):
        ind = cell_types == cell_type
        ax.scatter(xs[ind], ys[ind], zs[ind] + surface_z_offsets[zi],
                   edgecolors='none', marker='.', alpha=scatter_alpha,
                   s=marker_sz, zorder=zi,
                   color=cell_type_colors[cell_type], label=cell_type)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * box_width_scale, box.height])
    lgnd = ax.legend(loc='center left', fontsize=legend_font_sz, bbox_to_anchor=(box_to_anchor_x, box_to_anchor_y),
                     scatterpoints=1, handletextpad=0.0,
                     handlelength=1.25, borderaxespad=.0,
                     alignment='left', columnspacing=.5)

    for handle in lgnd.legendHandles:
        handle._sizes = [legend_handle_sz]

    ax.plot_surface(xx, yy, zz + surface_z_offsets[zi], color=surface_color, alpha=surface_alpha, zorder=zi)
    ax.text(y_label_horizontal_offset, y_label_vertical_offset, ylabel_z_offsets[zi] + surface_z_offsets[zi], ylabels[zi], 'y',
            color=title_color, fontsize=title_font_sz, zorder=1e5, ha='left', va='center')

    zi = 2
    ccc_fp = f"data/{args.dataset}/{args.ccc_fn}"
    if not os.path.exists(ccc_fp):
        print(f"The CCC input file: {ccc_fp}, not exist, please double check the `ccc_fn` setting in config.py or your input")
        return -1
    df_net = pd.read_csv(ccc_fp, header=0, index_col=0)
    cts = df_net.index.to_numpy()
    cts = cts.astype(str) if cell_type_in_str else cts.astype(int)

    cell_type_coord_mapping = {ct: [x2s[cid], y2s[cid] + ccc_scatter_y_off] for cid, ct in enumerate(cts)}
    c2s = np.array([cell_type_colors[ct] for ct in cts])
    z2s = np.ones_like(x2s) * surface_z_offsets[zi]

    ax.scatter(x2s, y2s + ccc_scatter_y_off, z2s, c=c2s, s=ccc_scatter_sz, edgecolors='none', marker='.', alpha=ccc_scatter_alpha,
               zorder=zi)

    # add within-layer edges
    froms, tos = df_net.values.nonzero()
    colors = np.array([cell_type_colors[cts[fm]] for fm in froms])
    starts = np.array([cell_type_coord_mapping[cts[fm]] for fm in froms])
    ends = np.array([cell_type_coord_mapping[cts[to]] for to in tos])
    dists = np.array([np.linalg.norm(start - ends[sid]) for sid, start in enumerate(starts)])
    dists = dists / np.max(dists)

    np.random.seed(random_seed)

    for sid, start in enumerate(starts):
        if dists[sid] > 0.:
            exception_condition = arrow_exception and cts[froms[sid]] == 4 and cts[tos[sid]] == 6
            d = np.linalg.norm(start - ends[sid]) * (0.07 if exception_condition else ccc_d)
            direction = 1 if exception_condition else None
            x, y = draw_curve(start, ends[sid], t=ccc_t, d=d, random_arrow_direction=random_arrow_direction, direction=direction)
            c = arrow_path(x, y, t=ccc_t, push=ccc_push)
            ax.plot(x, y, [surface_z_offsets[zi]] * len(x), color=colors[sid], linewidth=lw, alpha=1)  # linewidths[sid]
            for d, dt in c:
                ax.arrow3D(d[0], d[1], surface_z_offsets[zi],
                           dt[0], dt[1], 0, fc=colors[sid], ec=colors[sid],
                           mutation_scale=mutation_scale)  # linewidths[sid]*3.
        else:
            centX = start[0] + ccc_center_x
            centY = start[1] + ccc_center_y
            theta = np.linspace(0, 2 * np.pi, 201)
            radius = ccc_radius
            x = centX + radius * np.cos(theta)
            y = centY + radius * np.sin(theta)
            ax.plot(x, y, [surface_z_offsets[zi]] * len(theta), color=colors[sid], linewidth=lw, alpha=1)
            c = arrow_path(x, y, t=ccc_t, push=ccc_self_loop_push)

            for d, dt in c:
                ax.arrow3D(d[0], d[1], surface_z_offsets[zi],
                           dt[0], dt[1], 0, fc=colors[sid], ec=colors[sid],
                           mutation_scale=mutation_scale)

    ax.plot_surface(xx, yy, zz + surface_z_offsets[zi], color=surface_color, alpha=surface_alpha, zorder=zi)
    # add label

    ax.text(y_label_horizontal_offset, y_label_vertical_offset, ylabel_z_offsets[zi] + surface_z_offsets[zi], ylabels[zi], 'y',
            color=title_color, fontsize=title_font_sz, zorder=1e5, ha='left', va='center')

    # set them all at the same x,y,zlims
    ax.set_ylim(min(ys) - ydiff * 0.05, max(ys) + ydiff * 0.05)
    ax.set_xlim(min(xs) - xdiff * xlim_l_ext_offset, max(xs) + xdiff * xlim_r_ext_offset)
    ax.set_zlim(zlim[0], zlim[1])

    ax.view_init(view_height, view_angle)
    ax.dist = zoom_level
    ax.set_axis_off()
    if return_ax:
        return ax
    else:
        mkdir(os.path.dirname(args.fig_save_fp))
        plt.savefig(args.fig_save_fp, dpi=425, bbox_inches='tight')