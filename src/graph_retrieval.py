import os
import subprocess
from io import StringIO
import pandas as pd
# import pygraphviz
import networkx as nx
from networkx.drawing.nx_agraph import read_dot
import matplotlib.pyplot as plt


def get_crispr_df(fsource):
    cmd = "crisprtools stat -aH " + fsource
    p = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
    fbuff = StringIO(p.stdout.decode('utf-8'))
    df = pd.read_csv(fbuff, sep="\t")
    return df


def get_crispr_df_4fun(v):
    df = pd.read_csv('../data/test/crispr.df.csv', index_col=0)
    q = df[df.GID == v]
    seq = None
    if 0 < q.shape[0] < 2:
        seq = q.iloc[0, 1]
    return seq

def draw_graph_from_dot(dot_fpath):
    G = read_dot(dot_fpath)
    plt.figure(figsize=(10, 10), facecolor="w", frameon=False)
    pos = nx.spring_layout(G)
    offset = .05
    pos_labels = {}
    keys = pos.keys()
    for key in keys:
        x, y = pos[key]
        pos_labels[key] = (x, y + offset)

    node_colors = [G.node[key]['fillcolor'] for key in G.node]
    nx.draw_networkx_nodes(G, pos, node_size=1200, node_shape='o', node_color=node_colors)
    nx.draw_networkx_edges(G, pos, width=2, edge_color='b')
    nx.draw_networkx_labels(G, pos_labels, fontsize=2)
    return


def plot_spacer_graph(dpath, gid, seq):
    dot_fpath = dpath + "Spacers_" + gid[1:] + "_" + seq + "_spacers.gv"
    draw_graph(dot_fpath)
    plt.show()
    return


def plot_spacer_graph_by_ID(crispr_path, graph_path, gid):
    """
    Plots the spacer graph on a Jupyter notebook.
    :param crispr_path: Path to the .crispr file.
    :param graph_path: Path to the directory containing the spacer graphs.
    :param gid: The group id of the spacer graph.
    :return:
    """
    df = get_crispr_df(crispr_path)
    q = df[df.GID == gid]
    seq = None
    if 0 < q.shape[0] < 2:
        seq = q.iloc[0, 1]
        plot_spacer_graph(graph_path, gid, seq)
    else:
        print("The graph", gid, "could not be identified on the .crispr file!")
    return


def save_spacer_graph(dpath, gid, seq, dout):
    dot_fpath = dpath + "Spacers_" + gid[1:] + "_" + seq + "_spacers.gv"
    if os.path.isfile(dot_fpath):
        draw_graph(dot_fpath)
        plt.savefig(dout + gid + '.png')
        plt.close("all")  # clearing the RAM
        print("Graph", gid, "saved at", dout)
    else:
        print("File not found at", dot_fpath)
    return


def search_box():
    import ipywidgets as ipyw
    from IPython.display import display
    tb = ipyw.widgets.Text(
        description='GID:',
        value='G37640',
    )
    tbpath = ipyw.widgets.Text(
        description='Path:',
        value='../data/test/',
    )
    bt = ipyw.widgets.Button(description='Search')
    box = ipyw.HBox([tb, tbpath, bt])
    display(box)

    def on_button_clicked(b):
        print("Searching for", tb.value)
        seq = get_crispr_df_4fun(tb.value)
        if seq is not None:
            plot_spacer_graph(tbpath.value, tb.value, seq)
        else:
            print('No spacer graph found!')

    bt.on_click(on_button_clicked)
    return


def generate_all_graphs(fsource, dpath, dout):
    """
    Generate all spacer graph plots, given a .crispr file source.
    :param fsource: The .crispr file path.
    :param dpath: Path to the directory where the spacer .gv graphs are located.
    :param dout: Path to the directory where the graph images must be generated.
    :return:
    """
    if not os.path.exists(dout):
        os.mkdir(dout)
    df = get_crispr_df(fsource)
    gids = df.GID.values  # get the graph ids from the crispr file
    for gid in gids:
        q = df[df.GID == gid]
        seq = None
        if 0 < q.shape[0] < 2:
            seq = q.iloc[0, 1]
        save_spacer_graph(dpath, gid, seq, dout)
    return
