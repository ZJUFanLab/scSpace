from utils import *
from models import MLP, TCA
from config import loadArgums
import torch


def main():
    print("*************** scSpace *****************")
    args = loadArgums()

    # load data
    sc_meta, sc_data, st_meta, st_data = load_data(args.project_name)

    print("=== Transfer Component Analysis ===")
    tca = TCA.TCA(kernel_type=args.kernel_type, dim=args.dim, lamb=args.lamb, gamma=args.gamma)
    sc_new, st_new = tca.fit(Xs=np.array(sc_data.T), Xt=np.array(st_data.T))
    print("Transfer Component Analysis done.")

    # load data and model
    st_dataloader = get_data_loader(args, st_meta=st_meta, st_data_new=st_new)
    if args.sample_mlp:
        mlp = init_model(net=MLP.sample_MLPEncoder(input_size=args.dim, hidden_size=args.hidden_size))
    else:
        mlp = init_model(net=MLP.MLPEncoder(input_size=args.dim))

    print("=== Training encoder for source domain ===")
    losses = []
    encoder = train(args, encoder=mlp, data_loader=st_dataloader, losses=losses)

    predict = encoder(make_cuda(torch.tensor(sc_new).to(torch.float32))).cpu().detach().numpy()
    pd.DataFrame(predict).to_csv(osp.join(args.project_name, 'processed_data/pseudo_space.csv'))


if __name__ == "__main__":
    main()