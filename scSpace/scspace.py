import numpy as np

from .models import *
from .utils import *
from natsort import natsorted


def construct_pseudo_space(
        sc_adata,
        st_adata,
        kernel_type: str = 'primal',
        dim: int = 50,
        lamb: int = 1,
        gamma: int = 1,
        batch_size: int = 16,
        hidden_size: int = 128,
        common_size: int = 2,
        activation: str = 'sigmoid',
        lr: float = 0.001,
        epoch_num: int = 1000,
        log_epoch: int = 100
):
    sc_data = sc_adata.X
    st_data = st_adata.X
    print('Beginning Transfer Component Analysis...')
    tca = TCA(kernel_type=kernel_type, dim=dim, lamb=lamb, gamma=gamma)
    sc_new, st_new = tca.fit(Xs=sc_data, Xt=st_data)
    print("Transfer Component Analysis done.")
    sc_adata.obsm['TCA'] = sc_new
    st_adata.obsm['TCA'] = st_new

    st_dataloader = get_data_loader(st_adata=st_adata, batch_size=batch_size)
    mlp = init_model(net=sample_MLPEncoder(input_size=dim, common_size=common_size, hidden_size=hidden_size,
                                           activation=activation))
    # mlp = init_model(net=MLPEncoder(input_size=dim))
    print('Beginning training encoder for source domain...')
    losses = []
    encoder = train(encoder=mlp, losses=losses, data_loader=st_dataloader,
                    lr=lr, epoch_num=epoch_num, log_epoch=log_epoch)
    print('Encoder for source domain training finished.')

    encoder.eval()
    sc_tca = np.array(sc_adata.obsm['TCA'])
    predict = encoder(make_cuda(torch.tensor(sc_tca).to(torch.float32)))
    pseudo_space = predict.cpu().detach().numpy()
    sc_adata.obsm['pseudo_space'] = pseudo_space

    return sc_adata, st_adata


def spatial_cluster(
        sc_adata,
        use_neighbors: bool = True,
        l: float = 10,
        Ks: int = 10,
        Kg: int = 20,
        n_comps: int = 50,
        alpha: int = 0,
        beta: int = 0,
        res: float = 0.5,
        target_num: Optional[int] = None,
        step: float = 0.1,
        max_run: int = 20,
        random_seed: int = 123,
        use_weights: bool = True,
        partition_type: Optional[Type[MutableVertexPartition]] = None,
        n_iterations: int = -1
):

    # create space graph
    _, knn_idx = knn(data=sc_adata.obsm['pseudo_space'], query=sc_adata.obsm['pseudo_space'], k=Ks+1)
    knn_idx = np.delete(knn_idx, 0, axis=1)
    N = knn_idx.shape[0]
    W = np.diag([1] * N)
    for i in range(N):
        W[i, knn_idx[i]] = 1

    # create gene expression graph
    sc_adata_raw = sc_adata.raw.to_adata()
    sc_adata_raw = sc_adata_raw[:, sc_adata_raw.var.highly_variable]
    # sc.tl.pca(sc_adata_raw, svd_solver='arpack')
    sc.tl.pca(sc_adata_raw, svd_solver='arpack', n_comps=n_comps)
    sc_adata.obsm['X_pca'] = sc_adata_raw.obsm['X_pca']

    _, knn_idx_2 = knn(data=sc_adata.obsm['X_pca'], query=sc_adata.obsm['X_pca'], k=Kg+1)
    knn_idx_2 = np.delete(knn_idx_2, 0, axis=1)

    if use_neighbors:
        print('Calculating spatial weights using the neighbours of spots...')
        # create spatial-weighted gene expression graph
        sp_weight = np.zeros((knn_idx_2.shape[0], knn_idx_2.shape[1]))
        sp_graph = ig.Graph.Adjacency((W > 0).tolist())
        sp_graph = ig.Graph.simplify(sp_graph)
        for i in range(knn_idx_2.shape[0]):
            to_idx = knn_idx_2[i, ]
            sp_weight[i, ] = sp_graph.shortest_paths(i, to_idx, mode='all')[0]

        sp_weight = 1 / (alpha + sp_weight) + beta
        adjacency = np.zeros((sp_weight.shape[0], sp_weight.shape[0]))
        for i in range(N):
            adjacency[i, knn_idx_2[i]] = sp_weight[i]
    else:
        print('Calculating spatial weights using the distances across spots...')
        pairwise_dist = euclid_dist(coord=sc_adata.obsm['pseudo_space'])
        sp_weight = np.exp(-1 * (pairwise_dist**2) / (2 * (l**2)))
        adjacency = np.zeros((knn_idx_2.shape[0], knn_idx_2.shape[0]))
        for i in range(knn_idx_2.shape[0]):
            to_idx = knn_idx_2[i, ]
            adjacency[i, knn_idx_2[i]] = sp_weight[i][to_idx]

    g = get_igraph_from_adjacency(adjacency=adjacency, directed=None)

    if target_num is None:
        print('Unsupervised clustering with res = ', res, '...')
        cluster = run_leiden(g=g, resolution=res, random_state=random_seed, use_weights=use_weights,
                             n_iterations=n_iterations, partition_type=partition_type)
    else:
        assert target_num > 1, 'Target cluster number must be greater than 1!'
        recommended_res = search_res(g=g, target_num=target_num, start=res, step=step,
                                     random_state=random_seed, max_run=max_run)
        cluster = run_leiden(g=g, resolution=recommended_res, random_state=random_seed, use_weights=use_weights,
                             n_iterations=n_iterations, partition_type=partition_type)

    sc_adata.obs['scSpace'] = pd.Categorical(
        values=cluster.astype('U'),
        categories=natsorted(map(str, np.unique(cluster))),
    )

    return sc_adata


# def spatial_cluster_sub(
#         sc_adata,
#         celltype_key: str,
#         select_cell_type: list,
#         Ks: int = 10,
#         Kg: int = 20,
#         alpha: int = 0,
#         beta: int = 0,
#         res: float = 0.5,
#         target_num: Optional[int] = None,
#         step: float = 0.1,
#         random_seed: int = 123,
#         use_weights: bool = True,
#         partition_type: Optional[Type[MutableVertexPartition]] = None,
#         n_iterations: int = -1
# ):
#     sc_adata_sub = sc_adata[sc_adata.obs[celltype_key].isin(select_cell_type)]
#     sc_adata_sub = spatial_cluster(sc_adata=sc_adata_sub, Ks=Ks, Kg=Kg, alpha=alpha, beta=beta, res=res,
#                                    target_num=target_num, step=step, random_seed=random_seed, use_weights=use_weights,
#                                    partition_type=partition_type, n_iterations=n_iterations)
#
#     sc_adata.obs['scSpace'] = sc_adata.obs[celltype_key]
#     sc_adata.obs.loc[sc_adata_sub.obs.index, 'scSpace'] = sc_adata_sub.obs['scSpace']
#
#     return sc_adata


