import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import rankdata
from adjustText import adjust_text


def tf_dotplot(adata,
               tfs_to_plot=None,
               group="cell_type",
               color="tf",
               size="tf_sig"
               ):

    color_df=adata.obsm[color].groupby(adata.obs[group]).mean()
    size_df=adata.obsm[size].groupby(adata.obs[group]).mean()

    if tfs_to_plot is None:
        tfs_to_plot=adata.obsm[color].index

    tf_ct_mean=color_df.divide(
        np.sqrt(np.square(color_df).sum(axis=0)),
        axis=1)

    plot_df=pd.melt(tf_ct_mean.reset_index(), id_vars=group)

    plot_df['p_val']=pd.melt(
        size_df.reset_index(), id_vars=group)['value']
    plot_df.columns=["Cell_Type","TF",  "TF_Activity", "Proportion_Significant"]

    plt.figure(figsize=(3,0.33*len(tfs_to_plot)))

    ax =sns.scatterplot(
        data=plot_df.query("TF in @tfs_to_plot"),
        y="TF", x="Cell_Type", hue="TF_Activity", size="Proportion_Significant",
        palette="PiYG",sizes=(0, 250)
    )
    ax.set_ylim(-0.5, -0.5+len(tfs_to_plot))
    ax.set_xticklabels(size_df.index,rotation = 90)
    plt.legend(bbox_to_anchor=(1.05,1),
               loc='upper left',
               borderaxespad=0)
def tf_protien_line_plot(tf_protein, protein, title=None, cutoff=0.8 ):
    cors=tf_protein[protein]
    fig, ax = plt.subplots(figsize= (12,9))

    ax.scatter(
        rankdata(cors.to_list()),cors, s=50, alpha=0.8
    )

    TEXTS = []
    for i in range(len(cors)):

        if abs(cors.iloc[i])>cutoff:
            x = rankdata(cors.to_list())[i]
            y = cors.iloc[i]
            text = cors.index[i]
            TEXTS.append(ax.text(x, y, text, fontsize=16, fontname="Poppins"))


    # Adjust text position and add arrows ----------------------------
    # 'expand_points' is a tuple with two multipliers by which to expand
    # the bounding box of texts when repelling them from points

    # 'arrowprops' receives a dictionary with all the properties we want
    # for the arrows
    adjust_text(
        TEXTS,
        expand_points=(3,3),
        arrowprops=dict(
            arrowstyle="->",
            lw=2,
            color='b',
            alpha=0.5
        ),
        only_move={ 'text':'y', 'objects':'y'},
        force_text=2,
        ax=fig.axes[0]
    )
    ax.set_ylabel("Correlation", fontdict={"size": 16})
    ax.set_xlabel("Transcription Factor", fontdict={"size": 16})
    ax.set_title(title,fontdict={"size": 18})

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    return fig