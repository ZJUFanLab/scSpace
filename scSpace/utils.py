import math
from typing import Optional, Tuple, Sequence, Type
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from torch import nn
import torch
import torch.backends.cudnn as cudnn
import torch.optim as optim
from sklearn.neighbors import KDTree
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
import igraph as ig
import leidenalg
import logging as logg

try:
    from leidenalg.VertexPartition import MutableVertexPartition
except ImportError:

    class MutableVertexPartition:
        pass


    MutableVertexPartition.__module__ = 'leidenalg.VertexPartition'


class MyDataset(Dataset):  # 需要继承data.Dataset
    def __init__(self, data, label):
        # 1. Initialize file path or list of file names.
        self.data = data
        self.label = label

    def __getitem__(self, idx):
        # TODO
        tmp_x = self.data[idx]
        tmp_y_tag = self.label[idx]

        return (tmp_x, tmp_y_tag)  # tag 分类

    def __len__(self):
        # You should change 0 to the total size of your dataset.
        return self.label.shape[0]


def load_data(
        sc_data_path: str,
        sc_meta_path: str,
        st_data_path: str,
        st_meta_path: str,
        spatial_key: list = ['xcoord', 'ycoord']
):
    """
    load scRNA-seq data and spatial transcriptomic data, and return to AnnData objects
    :param sc_data_path: input scRNA-seq data file path
    :param sc_meta_path: input scRNA-seq meta file path
    :param st_data_path: input spatial transcriptomic data file path
    :param st_meta_path: input spatial transcriptomic data file path
    :param spatial_key: the columns of spatial coordinates
    :return: AnnData objects of scRNA-seq data and spatial transcriptomic data
    """

    print('Loading data...')
    sc_meta = pd.read_csv(sc_meta_path, index_col=0)
    sc_data = pd.read_csv(sc_data_path, index_col=0)
    st_meta = pd.read_csv(st_meta_path, index_col=0)
    st_data = pd.read_csv(st_data_path, index_col=0)
    print('Data have been loaded.')

    sc_adata = ad.AnnData(sc_data.T)
    sc_adata.obs = sc_meta

    st_adata = ad.AnnData(st_data.T)
    st_adata.obs = st_meta
    st_adata.obsm['spatial'] = np.array(st_meta[spatial_key])

    return sc_adata, st_adata


# sc_obj, st_obj = load_data(sc_data_path='/home/qjy321/workspace/scspace_dev/scSpace/data/demo_sc_data.csv',
#                            sc_meta_path='/home/qjy321/workspace/scspace_dev/scSpace/data/demo_sc_meta.csv',
#                            st_data_path='/home/qjy321/workspace/scspace_dev/scSpace/data/demo_st_data.csv',
#                            st_meta_path='/home/qjy321/workspace/scspace_dev/scSpace/data/demo_st_meta.csv')


def preporcess(
        sc_adata,
        st_adata,
        st_type: str = 'spot',
        n_features: int = 2000,
        normalize: bool = True,
        select_hvg: str = 'intersection'
):
    """
    pre-process the scRNA-seq and spatial transcriptomics data (find HVGs and normalized the data)
    :param select_hvg: 'intersection' or 'union'
    :param sc_adata: AnnData object of scRNA-seq data
    :param st_adata: AnnData object of spatial transcriptomics data
    :param st_type: the type of spatial transcriptomics data, `spot` or `image`
    :param n_features: the number of HVGs to select
    :param normalize: whether to normalize the data or not
    :return: AnnData object of processed scRNA-seq data and spatial transcriptomic data
    """

    assert sc_adata.shape[1] >= n_features, 'There are too few genes in scRNA-seq data, please check again!'
    sc.pp.highly_variable_genes(sc_adata, flavor="seurat_v3", n_top_genes=n_features)

    assert st_type in ['spot',
                       'image'], 'Please select the correct type of spatial transcriptomic data, `spot` or `image`!'

    if st_type == 'spot':
        assert st_adata.shape[1] >= n_features, 'There are too few genes in ST data, please check again!'
        sc.pp.highly_variable_genes(st_adata, flavor="seurat_v3", n_top_genes=n_features)
    elif st_type == 'image':
        if st_adata.shape[1] >= n_features:
            sc.pp.highly_variable_genes(st_adata, flavor="seurat_v3", n_top_genes=n_features)
        else:
            sc.pp.highly_variable_genes(st_adata, flavor="seurat_v3", n_top_genes=st_adata.shape[1])

    if normalize:
        # sc_adata
        sc.pp.normalize_total(sc_adata, target_sum=1e4)
        sc.pp.log1p(sc_adata)
        # st_adata
        sc.pp.normalize_total(st_adata, target_sum=1e4)
        sc.pp.log1p(st_adata)

    sc_adata.raw = sc_adata
    st_adata.raw = st_adata

    sc_hvg = sc_adata.var['highly_variable'][sc_adata.var['highly_variable'] == True].index
    st_hvg = st_adata.var['highly_variable'][st_adata.var['highly_variable'] == True].index
    if select_hvg == 'intersection':
        inter_gene = set(sc_hvg).intersection(set(st_hvg))
    elif select_hvg == 'union':
        sc_gene = set(sc_adata.var_names)
        st_gene = set(st_adata.var_names)
        common_gene = set(sc_gene).intersection(set(st_gene))

        inter_gene = set(sc_hvg).union(set(st_hvg))
        inter_gene = set(inter_gene).intersection(set(common_gene))

    sc_adata = sc_adata[:, list(inter_gene)]
    st_adata = st_adata[:, list(inter_gene)]

    print('Data have been pre-processed.')

    return sc_adata, st_adata


# sc_obj, st_obj = preporcess(sc_adata=sc_obj, st_adata=st_obj, st_type='image', n_features=2000, normalize=True)


def init_model(net):
    """Init models with cuda."""

    # check if cuda is available
    if torch.cuda.is_available():
        cudnn.benchmark = True
        net.cuda()

    return net


def make_cuda(tensor):
    """Use CUDA if it's available."""
    if torch.cuda.is_available():
        tensor = tensor.cuda()
    return tensor


def get_data_loader(
        st_adata,
        batch_size: int = 16,
):
    """Get dataset loader."""
    # dataset
    data = torch.tensor(np.array(st_adata.obsm['TCA'])).to(torch.float32)
    label = torch.tensor(np.array(st_adata.obsm['spatial'])).to(torch.float32)

    # data loader
    dataloader = DataLoader(
        dataset=MyDataset(data=data, label=label),
        batch_size=batch_size,
        shuffle=True)

    return dataloader


def train(
        encoder,
        data_loader,
        losses,
        lr: float = 0.001,
        epoch_num: int = 1000,
        log_epoch: int = 100
):
    """Train source domain."""
    ####################
    # 1. setup network #
    ####################

    # set train state for Dropout and BN layers
    encoder.train()

    # setup criterion and optimizer
    optimizer = optim.Adam(encoder.parameters(), lr=lr)
    criterion = nn.MSELoss(reduction='mean')

    ####################
    # 2. train network #
    ####################

    for epoch in tqdm(range(epoch_num)):
        batch_loss = []
        for step, data in enumerate(data_loader):
            features, labels = data

            features = make_cuda(features)
            labels = make_cuda(labels)

            # zero gradients for optimizer
            optimizer.zero_grad()

            # compute loss for critic
            pred = encoder(features)
            loss = criterion(pred, labels)

            # optimize source classifier
            loss.backward()
            optimizer.step()

            batch_loss.append(loss.cpu().data.numpy())

        # print
        if (epoch + 1) % log_epoch == 0:
            losses.append(np.mean(batch_loss))
            print("Epoch [{}/{}]: Batch loss={}"
                  .format(epoch + 1,
                          epoch_num,
                          np.mean(batch_loss)))

    return encoder


def knn(data, query, k):
    tree = KDTree(data)
    dist, ind = tree.query(query, k)
    return dist, ind


def euclid_dist(coord):
    """
    Calculate Euclidean distance between points
    """
    dis = []
    for i in range(len(coord)):
        dis.append([])
        for j in range(len(coord)):
            dis[i].append(math.dist(coord[i], coord[j]))
    dis = np.array(dis)
    return dis


def get_igraph_from_adjacency(adjacency, directed=None):
    """Get igraph graph from adjacency matrix."""

    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except KeyError:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warning(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g


def run_leiden(
        g,
        resolution: float = 0.5,
        random_state: int = 123,
        use_weights: bool = True,
        n_iterations: int = -1,
        partition_type: Optional[Type[MutableVertexPartition]] = None,
):
    """
    run Leiden clustering algorithm
    :param g: spatial-weighted gene expression graph
    :param resolution: A parameter value controlling the coarseness of the clustering.
                       Higher values lead to more clusters.
    :param random_state: Change the initialization of the optimization.
    :param use_weights: If `True`, edge weights from the graph are used in the computation
    :param n_iterations: How many iterations of the Leiden clustering algorithm to perform.
                         Positive values above 2 define the total number of iterations to perform,
                         -1 has the algorithm run until it reaches its optimal clustering.
    :param partition_type: Type of partition to use. Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
                           For the available options, consult the documentation for :func:`~leidenalg.find_partition`.
    :return: clustering results
    """

    partition_kwargs = dict()

    if partition_type is None:
        partition_type = leidenalg.RBConfigurationVertexPartition
    if use_weights:
        partition_kwargs['weights'] = np.array(g.es['weight']).astype(np.float64)
    partition_kwargs['n_iterations'] = n_iterations
    partition_kwargs['seed'] = random_state
    if resolution is not None:
        partition_kwargs['resolution_parameter'] = resolution
    # clustering proper
    part = leidenalg.find_partition(g, partition_type, **partition_kwargs)
    groups = np.array(part.membership)

    return groups


def search_res(
        g,
        target_num,
        start: float = 0.5,
        step: float = 0.1,
        random_state: int = 123,
        max_run: int = 10,

):
    res = start
    print("Start at res = ", res, "step = ", step)
    original_group = run_leiden(g=g, resolution=res, random_state=random_state)
    original_num = len(np.unique(original_group))
    print("Res = ", res, "number of clusters = ", original_num)
    run = 0

    while original_num != target_num:
        sign = 1 if (original_num < target_num) else -1
        new_group = run_leiden(g=g, resolution=(res + step * sign), random_state=random_state)
        new_num = len(np.unique(new_group))
        print("Res = ", (res + step * sign), "number of clusters = ", new_num)

        if new_num == target_num:
            res = res + step * sign
            print("Recommended res = ", res)
            return res

        new_sign = 1 if (new_num < target_num) else -1
        if new_sign == sign:
            res = res + step * sign
            print("Res changed to", res)
            original_num = new_num
        else:
            step = step / 2
            print("Step changed to", step)

        assert res > 0, 'Res must be positive! Please increase the number of target clusters!'

        if run > max_run:
            print('No resolution found, recommended res = ', res)
            return res
        run += 1

    print("Recommended res = ", res)

    return res


# TODO: Add comments
def cal_dist(coord,
             normalize: bool = True):
    dist = []
    for i in tqdm(range(len(coord))):
        xi, yi = coord[i, :]
        for j in range(len(coord)):
            if i >= j:
                continue
            xj, yj = coord[j, :]
            dist_tmp = np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
            dist.append(dist_tmp)

    if normalize:
        dist = dist / max(dist)

    return dist


def cal_dist_group(sc_adata, group_key, select_group):
    coord_all = sc_adata.obsm['pseudo_space']

    group_df = sc_adata.obs[group_key]
    group_df = group_df.reset_index()
    group = list(np.unique(group_df[group_key]))
    assert select_group in group, 'Please select the correct group!'

    select_idx = group_df[group_df[group_key] == select_group].index
    selsct_coord = coord_all[select_idx]

    #     dist_all = pd.DataFrame(columns=['dist', 'group'])
    dist_all = []
    for g in group:
        print('Calculating all cell pairs between', select_group, 'and', g, '...')
        dist_group = []
        tmp_idx = group_df[group_df[group_key] == g].index
        tmp_coord = coord_all[tmp_idx]

        if g == select_group:
            for i in tqdm(range(len(tmp_coord))):
                xi, yi = tmp_coord[i, :]
                for j in range(len(selsct_coord)):
                    if i >= j:
                        continue
                    xj, yj = selsct_coord[j, :]
                    dist_tmp = np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
                    dist_group.append(dist_tmp)
        else:
            for i in tqdm(range(len(tmp_coord))):
                xi, yi = tmp_coord[i, :]
                for j in range(len(selsct_coord)):
                    xj, yj = selsct_coord[j, :]
                    dist_tmp = np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
                    dist_group.append(dist_tmp)

        dist_all.append(dist_group)

    return dist_all
