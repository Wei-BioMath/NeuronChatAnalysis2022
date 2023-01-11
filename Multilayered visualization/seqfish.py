import os
import math
import util
import numpy as np
import pandas as pd
from config import parser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
setattr(Axes3D, 'arrow3D', util.arrow3D)

if __name__ == "__main__":
    args = parser.parse_args()
    args.dataset = "seqfish"
    args.ann_img_fp = f"data/seqfish/cortex_svz_location_fields.png"
    args.fig_save_fp = f'figures/seqfish.pdf'
    meta_fp = f"data/{args.dataset}/{args.meta_fn}"
    ccc_fp = f"data/{args.dataset}/{args.ccc_fn}"
    if not os.path.exists(meta_fp) or not os.path.exists(ccc_fp):
        print(
            f"meta file: {meta_fp} or CCC input file: {ccc_fp}, not exist, please double check the `meta_fn` and `ccc_fn` settings in config.py or your input")
        exit(-1)

    df = pd.read_csv(meta_fp)
    df_net = pd.read_csv(ccc_fp, header=0, index_col=0)

    img = plt.imread(args.ann_img_fp)
    img_w, img_h, _ = img.shape
    img = np.flipud(img[:img_w // 2, int(img_w * .02):int(img_w * .8), :])

    cts = df_net.index.to_numpy().astype(str)
    x2s = np.array([np.mean(df.loc[df['cell_type'] == ct, "sdimx"].to_numpy().astype(float)) for ct in cts], dtype=float)
    y2s = np.array([np.mean(df.loc[df['cell_type'] == ct, "sdimy"].to_numpy().astype(float)) for ct in cts], dtype=float)
    random_seed=99
    ccc_scatter_y_off=500

    ax = util.plot_3d_spatial_ccc(args,
                             img,
                             x2s,
                             y2s,
                             view_angle=-80,
                             view_height=11,
                             zoom_level=8,
                             marker_sz=20,
                             lw=1.,
                             title_font_sz=7,
                             annotation_img_stride=1,
                             surface_xl_ext=.04,
                             surface_xr_ext=.06,
                             box_to_anchor_x=.865,
                             scatter_alpha=.8,
                             box_width_scale=1.1,
                             y_label_horizontal_offset=-1800,
                             y_label_vertical_offset=-1500,
                             surface_z_offsets=(.0, .4, .75),
                             ylabel_z_offsets= (0.08, .1, .12),
                             ccc_center_x=-150,
                             ccc_center_y=300.,
                             ccc_radius=300,
                             ccc_push=1.8,
                             ccc_d= .3,
                             ccc_self_loop_push=1.6,
                             ccc_scatter_y_off=ccc_scatter_y_off,
                             mutation_scale=5,
                             random_arrow_direction=False,
                             random_seed=random_seed,
                             return_ax=True,
                             )

    xs = df["sdimx"].to_numpy().astype(float)
    ys = df["sdimy"].to_numpy().astype(float)
    cell_types = df["cell_type"].to_numpy().astype(str)
    cell_type_coord_mapping = {ct: [x2s[cid], y2s[cid] + ccc_scatter_y_off] for cid, ct in enumerate(cts)}

    x_offsets = np.array([2250., 3780.97, 5805.72, 7280.07, 7755.57])
    y_offsets = np.array([-2600., -2600., -2600., -2600., -1750.])
    z_offsets = np.array([0.34, 0.34, 0.34, 0.34, 0.24])

    top_x_offset = list(x2s) + [x2s[-1]]
    top_y_offset = list(y2s) + [y2s[-1] + 1500]
    top_x_offset[0] = top_x_offset[0] - 500

    lines3d_between = [([x_offsets[i], y_offsets[i], .2], [top_x_offset[i], top_y_offset[i], z_offsets[i]]) for i
                       in range(len(x_offsets))]
    between_lines = Line3DCollection(lines3d_between, zorder=1, color='yellow',
                                     alpha=0.8, linestyle='--', linewidth=.5)
    ax.add_collection3d(between_lines)

    np.random.seed(random_seed)

    for cid, ctype in enumerate(np.unique(cell_types)):
        up_node = cell_type_coord_mapping[ctype]
        bottom_node_indexs = np.squeeze(np.where(cell_types == ctype))
        distances = [math.dist(up_node, [xs[node_index], ys[node_index]]) for node_index in bottom_node_indexs]
        min_dist_indices = np.argsort(distances)
        start = 22
        n_node_selection = 5
        step = 1
        lines3d_between = [
            ([xs[bottom_node_indexs[node_index]], ys[bottom_node_indexs[node_index]], .4], up_node + [.75]) for
            node_index in min_dist_indices[start: start + n_node_selection * step: step]]
        between_lines = Line3DCollection(lines3d_between, zorder=2, color='.5',
                                         alpha=0.4, linestyle='--', linewidth=.6)
        ax.add_collection3d(between_lines)

    util.mkdir(os.path.dirname(args.fig_save_fp))
    plt.savefig(args.fig_save_fp, dpi=425, bbox_inches='tight')