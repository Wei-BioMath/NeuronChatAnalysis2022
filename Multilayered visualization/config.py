# -*- coding:utf-8 -*-
import argparse

parser = argparse.ArgumentParser(description='Plot 3D Spatial Cell-Cell Communication Diagram')

parser.add_argument('-dataset', type=str, default="seqfish",
                    help='The name of the dataset used for plotting')

parser.add_argument('-meta_fn', type=str, default="meta.csv",
                    help='The meta filename for the dataset. Default: `meta.csv`')

parser.add_argument('-ccc_fn', type=str, default="net_top10.csv",
                    help='The filename of the top-n cell-cell communication \n'
                         'probability matrix. Default: `net_top10.csv`')

parser.add_argument('-ann_img_fp', type=str, default="./path_to_annotation_image.png",
                    help='The file path of the annotation image. \n'
                         'Default: `./path_to_annotation_image.png`')

parser.add_argument('-colormap', type=str, default="Spectral",
                    help='The colormap of scatter plot and cell-cell communication\n'
                         'plot. Default: `Spectral`, see more options at: \n'
                         '`https://matplotlib.org/stable/tutorials/colors/colormaps.html`')

parser.add_argument('-fig_save_fp', type=str, default="./path_to_save_fig.pdf",
                    help='The file path for the figure output. \n'
                         'Default: `./path_to_save_fig.pdf`')

