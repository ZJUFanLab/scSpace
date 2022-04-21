import numpy as np
import os.path as osp
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from torch import nn
import torch
import torch.backends.cudnn as cudnn
import torch.optim as optim


def load_data(data_path):
    sc_meta_path = osp.join(data_path, 'processed_data/sc_meta_processed.csv')
    sc_data_path = osp.join(data_path, 'processed_data/sc_data_processed.csv')
    st_meta_path = osp.join(data_path, 'processed_data/st_meta_processed.csv')
    st_data_path = osp.join(data_path, 'processed_data/st_data_processed.csv')
    print("Loading processed data......")

    sc_meta = pd.read_csv(sc_meta_path, index_col=0)
    sc_data = pd.read_csv(sc_data_path, index_col=0)
    st_meta = pd.read_csv(st_meta_path, index_col=0)
    st_data = pd.read_csv(st_data_path, index_col=0)
    print("Load data ok")

    return sc_meta, sc_data, st_meta, st_data


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


def get_data_loader(args, st_meta, st_data_new):
    """Get dataset loader."""
    # dataset
    data = torch.tensor(st_data_new).to(torch.float32)
    label = torch.tensor(np.array(st_meta[['xcoord', 'ycoord']])).to(torch.float32)

    # data loader
    dataloader = DataLoader(
        dataset=MyDataset(data=data, label=label),
        batch_size=args.batch_size,
        shuffle=True)

    return dataloader


def train(args, encoder, data_loader, losses):
    """Train source domain."""
    ####################
    # 1. setup network #
    ####################

    # set train state for Dropout and BN layers
    encoder.train()

    # setup criterion and optimizer
    optimizer = optim.Adam(encoder.parameters(), lr=args.lr)
    criterion = nn.MSELoss(reduction='mean')

    ####################
    # 2. train network #
    ####################

    for epoch in range(args.epoch_num):
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
        if (epoch + 1) % args.log_epoch == 0:
            losses.append(np.mean(batch_loss))
            print("Epoch [{}/{}]: Batch loss={}"
                  .format(epoch + 1,
                          args.epoch_num,
                          np.mean(batch_loss)))

    return encoder
