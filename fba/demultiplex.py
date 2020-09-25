# demulitplex.py

import re
import sys
import numpy as np
import pandas as pd
import scipy.sparse
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
from pyclustering.cluster.kmedoids import kmedoids
from sklearn.preprocessing import StandardScaler
from pathlib import Path
from sklearn.manifold import TSNE
from .utils import get_logger


logger = get_logger(logger_name=__name__)


params = {'pdf.fonttype': 42,
          'mathtext.default': 'regular',
          'axes.axisbelow': True
          }
plt.rcParams.update(params)


def normalize_clr(x):
    """Computes the centered log ratio (clr) transform."""

    return np.log1p(x / (np.exp(sum(np.log1p(x[x > 0])) / len(x))))


def get_cell_identity(x,
                      barcodes):
    """."""

    if sum(x) == 1:
        cell_identity = np.array(object=barcodes[x == 1], dtype=object)
    elif sum(x) == 0:
        cell_identity = np.array(object='negative',  dtype=object)
    elif sum(x):
        cell_identity = np.array(object='multiplet', dtype=object)

    return cell_identity


def cluster(m, seed=42):
    """."""

    np.random.seed(seed=seed)

    m = m.apply(normalize_clr, axis=1)
    kmedoids_instance = kmedoids(
        data=m.T.values,
        initial_index_medoids=np.random.randint(
            low=0,
            high=m.shape[1],
            size=m.shape[0] + 1,
            dtype=np.int),
        tolerance=0.001,
        ccore=True)

    kmedoids_instance.process()

    clusters_kmedoids = kmedoids_instance.get_clusters()

    return clusters_kmedoids


def demultiplex(m, seed=42, q=99.5):
    """."""

    m_identity = np.zeros_like(m)

    c = cluster(m, seed=seed)

    m_feature_avg = pd.DataFrame(
        data=np.zeros(
            shape=[m.shape[0], len(c)],
            dtype=np.float, order='C'),
        index=m.index)

    for index, value in enumerate(c):
        m_feature_avg.iloc[:, index] = m.iloc[:, value].mean(axis=1).values

    for index, value in enumerate(m_feature_avg.index):

        cells_selected = [i for j, i in enumerate(c)
                          if j != np.argmax(m_feature_avg.iloc[index].values)]
        cells_selected = [j for i in cells_selected for j in i]

        cells_selected_counts = m.iloc[index, cells_selected].values
        cells_selected_counts = cells_selected_counts[
            cells_selected_counts < np.percentile(a=cells_selected_counts,
                                                  q=q)]

        mod_nbin = sm.NegativeBinomial(endog=cells_selected_counts,
                                       exog=np.ones(
                                           len(cells_selected_counts)),
                                       loglike_method='nb2')
        res_nbin = mod_nbin.fit(disp=False)

        mu, alpha = res_nbin.params
        mu = np.exp(mu)
        size = 1 / alpha * mu ** 0
        prob = size / (size + mu)
        count_cutoff = scipy.stats.nbinom.ppf(q=0.99, n=size, p=prob)
        # print(value, ':', count_cutoff)

        m_identity[index, m.iloc[index, :] >= count_cutoff] = 1

        cells_demultiplexed = pd.Series(data=np.apply_along_axis(
            func1d=get_cell_identity,
            axis=0,
            arr=m_identity,
            barcodes=m.index.values).flatten(),
            index=m.columns
        ).to_frame(name='category')

    return cells_demultiplexed, m_identity


def prepare_heatmap_matrix(m_identity,
                           m_norm,
                           percentile_low=1,
                           percentile_high=99):
    """Prepares heatmap matrix for visualization."""

    num_positive = m_identity.sum(axis=0)
    cells_singlet = num_positive[num_positive == 1].index
    cells_multiplet = num_positive[num_positive > 1].index
    cells_negative = num_positive[num_positive == 0].index

    cells_singlet = m_identity.loc[:, cells_singlet].sort_values(
        by=list(m_identity.index),
        axis=1,
        ascending=False).columns.values

    cells_multiplet = m_identity.loc[:, cells_multiplet].sort_values(
        by=list(m_identity.index),
        axis=1,
        ascending=False).columns.values

    cells_negative = m_identity.loc[:, cells_negative].sort_values(
        by=list(m_identity.index),
        axis=1,
        ascending=False).columns.values

    heatmap_matrix = m_norm.loc[:, np.concatenate(
        (cells_singlet,
         cells_multiplet,
         cells_negative),
        axis=0)]

    heatmap_matrix = StandardScaler(
        copy=True,
        with_mean=True,
        with_std=True).fit_transform(heatmap_matrix.T)

    heatmap_limits = np.percentile(
        a=heatmap_matrix,
        q=[percentile_low,
           percentile_high])

    heatmap_matrix[heatmap_matrix < heatmap_limits[0]] = heatmap_limits[0]
    heatmap_matrix[heatmap_matrix > heatmap_limits[1]] = heatmap_limits[1]

    heatmap_matrix = pd.DataFrame(data=heatmap_matrix.T,
                                  # index=m_identity.index,
                                  # columns=m_identity.columns)
                                  index=m_norm.index,
                                  columns=m_norm.columns)

    return heatmap_matrix


def plot_heatmap_features_selected(heatmat_matrix,
                                   ax,
                                   color_map='plasma',
                                   title=None):
    """Plots heatmap of selected features across single cells."""

    sns.heatmap(data=heatmat_matrix,
                cmap=color_map,
                cbar=True,
                xticklabels=False,
                yticklabels=heatmat_matrix.index,
                cbar_kws={'orientation': 'vertical',
                          'pad': 0.025,
                          'shrink': 0.75,
                          'aspect': 15,
                          'label': 'Z score'},
                ax=ax)

    ax.tick_params(axis='both',
                   which='major',
                   direction='out',
                   width=0.8,
                   labelsize=7,
                   labelcolor='black',
                   colors='#333333')

    if title:
        ax.set_title(label=title, fontsize=8)

    # customize color bar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(axis='both',
                        which='major',
                        direction='out',
                        width=0.8,
                        labelsize=7,
                        labelcolor='black',
                        colors='#333333')

    cbar.set_label(label='Z score',
                   fontsize=8,
                   rotation=-90,
                   labelpad=10)

    return ax


def prepare_embeddings(cells,
                       m,
                       seed):
    """Embeds cells."""

    # t-SNE
    embedding_tsne = TSNE(n_components=2,
                          perplexity=30.0,
                          early_exaggeration=12.0,
                          learning_rate=200.0,
                          n_iter=3000,
                          n_iter_without_progress=300,
                          min_grad_norm=1e-07,
                          metric='euclidean',
                          init='random',
                          verbose=0,
                          random_state=seed,
                          method='barnes_hut',
                          angle=0.5).fit_transform(m.loc[:, cells].T)

    embedding = pd.DataFrame(embedding_tsne,
                             columns=['x_tsne', 'y_tsne'],
                             index=cells)

    return embedding


def plot_embedding(embedding,
                   ax,
                   title=None,
                   group=None,
                   show_group_labels=False,
                   marker_size=10):
    """Plots embedding."""

    if type(embedding) is not pd.DataFrame:
        embedding = pd.DataFrame(data=embedding[:, range(2)],
                                 columns=['x', 'y'])
    else:
        embedding = embedding.iloc[:, range(2)].copy()
        embedding.columns = ['x', 'y']

    if group is not None:
        if embedding.shape[0] == len(group):
            embedding['category'] = list(group)

            embedding['color'] = embedding['category'].map(
                dict(zip(sorted(embedding['category'].unique()),
                         sns.color_palette(
                             palette='husl',
                             n_colors=embedding['category'].nunique()))))

        else:
            embedding['category'] = np.array(np.isin(embedding.index,
                                                     group),
                                             dtype=np.int)

            embedding.sort_values(by=['category'], inplace=True)
            embedding['color'] = embedding['category'].map(
                dict(zip([0, 1], ['darkgrey', '#ff8c69'])))

    else:
        embedding['color'] = 'tab:blue'

    ax.scatter(x=embedding['x'],
               y=embedding['y'],
               s=marker_size,
               marker='.',
               c=embedding['color'],
               alpha=1,
               linewidths=0,
               rasterized=False,
               edgecolors=None)

    if group is not None and show_group_labels:
        group_labels = embedding.groupby('category').median()
        if type(show_group_labels) is not bool:
            group_labels = group_labels.loc[show_group_labels].copy()

        for i in group_labels.index:
            ax.annotate(text=str(i),
                        # s=str(i),
                        xy=(group_labels.loc[i]['x'],
                            group_labels.loc[i]['y']),
                        fontsize=7,
                        horizontalalignment='center',
                        verticalalignment='center')

    if title:
        ax.set_title(label=title, fontdict=None, loc='center', fontsize=8)

    for i in ['top', 'bottom', 'left', 'right']:
        ax.spines[i].set_linewidth(w=0.5)
        ax.spines[i].set_color(c='grey')

    # ax.set_axis_off()
    ax.xaxis.set_ticks(ticks=[])
    ax.yaxis.set_ticks(ticks=[])
    # ax.set_facecolor(color='white')

    return ax


def demultiplex_feature_barcoding(matrix_featurecount_file,
                                  output_directory='demultiplexed',
                                  seed=42):
    """."""

    output_directory = Path(output_directory)
    output_directory.mkdir(exist_ok=True)

    CELLS_DEMULTIPLEXED_FILE = output_directory / 'cells_demultiplexed.csv'
    CELLS_DEMULTIPLEXED_HEATMAP_PLOT = output_directory / \
        'Pyplot_heatmap_cells_demultiplexed.pdf'

    CELLS_DEMULTIPLEXED_EMBEDDING_FILE = output_directory / \
        'cells_embedding_demultiplexed.csv'
    CELLS_DEMULTIPLEXED_EMBEDDING_PLOT = output_directory / \
        'Pyplot_embedding_cells_demultiplexed.pdf'

    logger.info(f'Output directory: {output_directory}')
    logger.info(f'Loading feature count matrix: {matrix_featurecount_file}')

    matrix_featurecount = pd.read_csv(
        filepath_or_buffer='matrix_featurecount.csv.gz',
        index_col=0
    )
    matrix_featurecount.index = [re.sub(
        pattern='_[A-Za-z]{1,}$', repl='', string=i)
        for i in matrix_featurecount.index]

    # matrix_featurecount = matrix_featurecount.loc[
    #     matrix_featurecount.sum(axis=1) >= 300]

    logger.info(f'Number of cells: {matrix_featurecount.shape[1]:,}')
    logger.info(f'Number of features: {matrix_featurecount.shape[0]:,}')
    logger.info(f'Total UMIs: {matrix_featurecount.values.sum():,}')
    logger.info('Median number of UMIs per cell: '
                + f'{np.median(matrix_featurecount.sum(axis=0)):,}')
    logger.info('Demultiplexing ...')

    try:
        cells_demultiplexed, m_identity = demultiplex(m=matrix_featurecount,
                                                      seed=seed,
                                                      q=99)
        cells_demultiplexed.to_csv(CELLS_DEMULTIPLEXED_FILE)
    except ValueError as err:
        logger.critical(
            'This demultiplexing method '
            + f'may not be the best solution for this dataset: {err}'
        )
        sys.exit(1)

    m_identity = pd.DataFrame(
        m_identity,
        index=matrix_featurecount.index,
        columns=matrix_featurecount.columns
    )

    # heatmap
    logger.info('Generating heatmap ...')
    matrix_heatmap = prepare_heatmap_matrix(
        m_identity=m_identity,
        m_norm=matrix_featurecount.apply(normalize_clr, axis=1),
        percentile_low=1,
        percentile_high=99
    )

    fig, ax = plt.subplots(
        nrows=1, ncols=1,
        figsize=(6, max(2, len(matrix_featurecount.index) * 0.25))
    )
    plot_heatmap_features_selected(heatmat_matrix=matrix_heatmap,
                                   ax=ax,
                                   color_map='viridis',
                                   title='Cell classification, kmedoids')
    plt.tight_layout()
    fig.savefig(fname=CELLS_DEMULTIPLEXED_HEATMAP_PLOT,
                transparent=True,
                bbox_inches='tight')

    # t-SNE
    logger.info('Embedding ...')
    cells_embedding = cells_demultiplexed.index[
        cells_demultiplexed['category'] != 'negative']

    embedding = prepare_embeddings(
        cells=cells_embedding,
        m=matrix_featurecount.apply(normalize_clr, axis=1),
        seed=seed)
    embedding['category'] = cells_demultiplexed.loc[embedding.index]
    embedding.to_csv(CELLS_DEMULTIPLEXED_EMBEDDING_FILE)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    plot_embedding(embedding=embedding,
                   ax=ax,
                   title='t-SNE; Cell classification, kmedoids',
                   group=embedding['category'],
                   show_group_labels=[
                       i for i in embedding['category'].unique()
                       if i != 'multiplet'],
                   marker_size=10)
    plt.tight_layout()

    fig.savefig(fname=CELLS_DEMULTIPLEXED_EMBEDDING_PLOT,
                transparent=None,
                bbox_inches='tight')

    return matrix_featurecount_file