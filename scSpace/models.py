import numpy as np
import scipy.io
import scipy.linalg
import sklearn.metrics
from torch import nn


def kernel(ker, X1, X2, gamma):
    K = None
    if not ker or ker == 'primal':
        K = X1
    elif ker == 'linear':
        if X2 is not None:
            K = sklearn.metrics.pairwise.linear_kernel(
                np.asarray(X1).T, np.asarray(X2).T)
        else:
            K = sklearn.metrics.pairwise.linear_kernel(np.asarray(X1).T)
    elif ker == 'rbf':
        if X2 is not None:
            K = sklearn.metrics.pairwise.rbf_kernel(
                np.asarray(X1).T, np.asarray(X2).T, gamma)
        else:
            K = sklearn.metrics.pairwise.rbf_kernel(
                np.asarray(X1).T, None, gamma)
    return K


class TCA:
    def __init__(self, kernel_type='primal', dim=30, lamb=1, gamma=1):
        '''
        Init func
        :param kernel_type: kernel, values: 'primal' | 'linear' | 'rbf'
        :param dim: dimension after transfer
        :param lamb: lambda value in equation
        :param gamma: kernel bandwidth for rbf kernel
        '''
        self.kernel_type = kernel_type
        self.dim = dim
        self.lamb = lamb
        self.gamma = gamma

    def fit(self, Xs, Xt):
        '''
        Transform Xs and Xt
        :param Xs: ns * n_feature, source feature
        :param Xt: nt * n_feature, target feature
        :return: Xs_new and Xt_new after TCA
        '''
        X = np.hstack((Xs.T, Xt.T))
        X /= np.linalg.norm(X, axis=0)
        m, n = X.shape
        ns, nt = len(Xs), len(Xt)
        e = np.vstack((1 / ns * np.ones((ns, 1)), -1 / nt * np.ones((nt, 1))))
        M = e * e.T
        M = M / np.linalg.norm(M, 'fro')
        H = np.eye(n) - 1 / n * np.ones((n, n))
        K = kernel(self.kernel_type, X, None, gamma=self.gamma)
        n_eye = m if self.kernel_type == 'primal' else n
        a, b = K @ M @ K.T + self.lamb * np.eye(n_eye), K @ H @ K.T
        w, V = scipy.linalg.eig(a, b)
        ind = np.argsort(w)
        A = V[:, ind[:self.dim]]
        Z = A.T @ K
        Z /= np.linalg.norm(Z, axis=0)

        Xs_new, Xt_new = Z[:, :ns].T, Z[:, ns:].T
        return Xs_new, Xt_new


activation_dict = {
    'relu': nn.ReLU,
    'sigmoid': nn.Sigmoid
}


class sample_MLPEncoder(nn.Module):
    """sample MLPEncoder encoder model for scSpace."""

    def __init__(self, input_size, common_size, hidden_size, activation):
        """Init MLP encoder."""
        super(sample_MLPEncoder, self).__init__()
        self.restored = False
        common_size = common_size
        self.linear = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            activation_dict[activation](),
            nn.Linear(hidden_size, common_size)
        )

    def forward(self, x):
        out = self.linear(x)
        return out


class MLPEncoder(nn.Module):
    """MLPEncoder encoder model for scSpace."""

    def __init__(self, input_size):
        """Init MLP encoder."""
        super(MLPEncoder, self).__init__()
        self.restored = False
        common_size = 2
        self.linear = nn.Sequential(
            nn.Linear(input_size, input_size // 2),
            nn.ReLU(inplace=True),
            nn.Linear(input_size // 2, input_size // 4),
            nn.ReLU(inplace=True),
            nn.Linear(input_size // 4, common_size)
        )

    def forward(self, x):
        out = self.linear(x)
        return out

