import os
import util
import numpy as np
import pandas as pd
from config import parser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
setattr(Axes3D, 'arrow3D', util.arrow3D)

if __name__ == "__main__":
    args = parser.parse_args()
    args.dataset = "visium"
    args.ann_img_fp = f"data/visium/mouse_brain_highres.png"
    args.fig_save_fp = f'figures/visium.pdf'

    meta_fp = f"data/{args.dataset}/{args.meta_fn}"
    ccc_fp = f"data/{args.dataset}/{args.ccc_fn}"

    if not os.path.exists(meta_fp) or not os.path.exists(ccc_fp):
        print(
            f"meta file: {meta_fp} or CCC input file: {ccc_fp}, not exist, please double check the `meta_fn` and `ccc_fn` settings in config.py or your input")
        exit(-1)

    df = pd.read_csv(meta_fp)
    df_net = pd.read_csv(ccc_fp, header=0, index_col=0)

    cts = df_net.index.to_numpy().astype(int)

    xys = np.array(
        [df.loc[df['cell_type'] == ct, ["sdimx", "sdimy"]].sample(random_state=47).to_numpy().astype(float).squeeze()
         for ct in cts], dtype=float)
    x2s, y2s = xys[:, 0], xys[:, 1]

    img = plt.imread(args.ann_img_fp)
    img_w, img_h, _ = img.shape
    img = np.flipud(img)

    util.plot_3d_spatial_ccc(args,
                             img,
                             x2s,
                             y2s,
                             view_angle=-75,
                             view_height=14,
                             zoom_level=8,
                             marker_sz=8,
                             title_font_sz=7,
                             annotation_img_stride=2,
                             surface_xl_ext=.05,
                             surface_xr_ext=.05,
                             box_to_anchor_x=.91,
                             scatter_alpha=.8,
                             box_width_scale=1.1,
                             y_label_horizontal_offset=-1200,
                             y_label_vertical_offset=-8000,
                             surface_z_offsets=(.0, .4, .75),
                             ylabel_z_offsets= (0.03, .06, .09),
                             ccc_center_x=100,
                             ccc_center_y=-500.,
                             ccc_radius=500,
                             ccc_push=1.,
                             ccc_self_loop_push=.9,
                             ccc_d=.03,
                             mutation_scale=5,
                             random_arrow_direction=True,
                             arrow_exception= True,
                             cell_type_in_str=False,
                             random_seed=25
                             )
