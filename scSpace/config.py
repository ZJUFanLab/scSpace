import argparse
import os.path as osp


def loadArgums():
    parser = argparse.ArgumentParser(description='scSpace')

    # project
    parser.add_argument('--project_name', '-pn', default='demo', help='project name', type=str)

    # Transfer Component Analysis
    parser.add_argument('--kernel_type', '-kt', default='primal', help='kernel type: primal/linear/rbf', type=str)
    parser.add_argument('--dim', '-d', default=50, help='dimension after transfer', type=int)
    parser.add_argument('--lamb', '-la', default=1, help='lambda value in equation', type=int)
    parser.add_argument('--gamma', '-ga', default=1, help='kernel bandwidth for rbf kernel', type=int)

    # Multilayer Perceptron
    parser.add_argument('--batch_size', '-bs', default=16, help='batch size eg:16', type=int)
    parser.add_argument('--lr', '-l', default=0.001, help='learning rate eg:0.001', type=float)
    parser.add_argument('--epoch_num', '-ep', default=1000, help='epoch number eg:1000', type=int)
    parser.add_argument('--log_epoch', '-le', default=100, help='log epoch eg:100', type=int)
    parser.add_argument("--sample_mlp", action='store_const', default=False, const=True)
    parser.add_argument('--hidden_size', '-hs', default=128, help='hidden_size for sample MLPEncoder eg:128', type=int)

    args = parser.parse_args()

    # path deal
    args.project_name = osp.abspath(osp.join(this_dir, 'data', args.project_name))

    return args


# path
this_dir = osp.dirname(__file__)
