from torch import nn


class sample_MLPEncoder(nn.Module):
    """sample MLPEncoder encoder model for scSpace."""

    def __init__(self, input_size, hidden_size):
        """Init MLP encoder."""
        super(sample_MLPEncoder, self).__init__()
        self.restored = False
        common_size = 2
        self.linear = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.Sigmoid(),
            nn.Linear(hidden_size, common_size)
        )

    def forward(self, x):
        out = self.linear(x)
        return out
