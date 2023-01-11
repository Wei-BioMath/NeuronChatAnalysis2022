import os
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
    args.dataset = "MERFISH"
    args.ann_img_fp = f"data/MERFISH/44_1.png"
    args.fig_save_fp = f'figures/MERFISH.pdf'
    meta_fp = f"data/{args.dataset}/{args.meta_fn}"
    ccc_fp = f"data/{args.dataset}/{args.ccc_fn}"
    if not os.path.exists(meta_fp) or not os.path.exists(ccc_fp):
        print(
            f"meta file: {meta_fp} or CCC input file: {ccc_fp}, not exist, please double check the `meta_fn` and `ccc_fn` settings in config.py or your input")
        exit(-1)

    df = pd.read_csv(meta_fp)
    df_net = pd.read_csv(ccc_fp, header=0, index_col=0)

    xs = df["sdimx"].to_numpy().astype(float)
    ys = df["sdimy"].to_numpy().astype(float)
    cts = df_net.index.to_numpy().astype(str)
    dx = xs.max() - xs.min()
    dy = ys.max() - ys.min()
    x2s = np.linspace(xs.min(), xs.max(), num=cts.shape[0])
    y2s = np.linspace(ys.min() + dx * .02, ys.max() - dx * .02, num=cts.shape[0])[::-1]

    img = plt.imread(args.ann_img_fp)
    img_w, img_h, _ = img.shape
    img = np.flipud(img)

    ax = util.plot_3d_spatial_ccc(args,
                             img,
                             x2s,
                             y2s,
                             view_angle=-80,
                             view_height=11,
                             zoom_level=8,
                             marker_sz=2,
                             title_font_sz=7,
                             annotation_img_stride=1,
                             surface_xl_ext=.05,
                             surface_xr_ext=.05,
                             box_to_anchor_x=.89,
                             box_to_anchor_y=.41,
                             scatter_alpha=.8,
                             box_width_scale=1.1,
                             y_label_horizontal_offset=-8300,
                             y_label_vertical_offset=-4900,
                             surface_z_offsets=(.0, .4, .75),
                             ylabel_z_offsets= (0., .03, .07),
                             ccc_center_x=100,
                             ccc_center_y=200.,
                             ccc_radius=180,
                             ccc_push=.8,
                             ccc_d= .2,
                             ccc_self_loop_push=.9,
                             random_arrow_direction=True,
                             random_seed=25,
                             return_ax=True
                             )

    xb = np.array([-4050, -4260, -3760, -4240])
    yb = np.array([-5650, -5200, -5320, -5600])
    xt = np.array([-3850, -6600, -3900, -6300])
    yt = np.array([-5650, -4200, -4200, -5300])
    zt = np.array([.4, .4, .45, .45])

    lines3d_between = [([xb[i], yb[i], .24], [xt[i], yt[i], zt[i]]) for i in range(len(xb))]
    between_lines = Line3DCollection(lines3d_between, zorder=1, color='.5',
                                     alpha=0.4, linestyle='--', linewidth=.5)
    ax.add_collection3d(between_lines)

    util.mkdir(os.path.dirname(args.fig_save_fp))
    plt.savefig(args.fig_save_fp, dpi=425, bbox_inches='tight')