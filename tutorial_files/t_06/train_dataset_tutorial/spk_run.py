#!/mms/common/conda/2022.05/envs/spk_ngn/bin/python3.9
import os
import torch
import logging
import schnetpack as spk
from schnetpack.utils import (
    get_dataset,
    get_metrics,
    get_loaders,
    get_statistics,
    get_model,
    get_trainer,
    ScriptError,
    evaluate,
    setup_run,
    get_divide_by_atoms,
#    get_indices_mlmm,
)
from schnetpack.utils.script_utils.settings import get_environment_provider
from schnetpack.utils.script_utils.parsing import build_parser

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
def run_prediction(model, loader, device, args):
    from tqdm import tqdm
    import numpy as np

    predicted = {}
    qm_values = {}
    for batch in tqdm(loader, ncols=120):
        batch = {
            k: v.to(device)
            for k, v in batch.items()
        }
        result = model(batch)
        for prop in result:
            if prop in predicted:
                predicted[prop] += [result[prop].cpu().detach().numpy()]
            else:
                predicted[prop] = [result[prop].cpu().detach().numpy()]
        for prop in batch:
            if prop in qm_values:
                qm_values[prop] += [batch[prop].cpu().detach().numpy()]
            else:
                qm_values[prop] = [batch[prop].cpu().detach().numpy()]

    import numpy as np
    np.savez("prediction.npz",predicted)
    np.savez("reference.npz",qm_values)
 
    logging.info('Stored model predictions')



def evaluate_(args, model, train_loader, val_loader, test_loader, device, train_args):
    # Get property names from model
    properties=train_args.property
    header = ['Subset']
    metrics = []
    for prop in properties:
        header += [f'{prop}_MAE', f'{prop}_RMSE']
        metrics += [
            spk.metrics.MeanAbsoluteError(prop, prop),
            spk.metrics.RootMeanSquaredError(prop, prop)
        ]


    results = []
    if ('train' in args.split) or ('all' in args.split):
        logging.info('Training split...')
        results.append(['training'] + ['%.7f' % i for i in evaluate_dataset(metrics, model, train_loader, device,properties)])

    if ('validation' in args.split) or ('all' in args.split):
        logging.info('Validation split...')
        results.append(['validation'] + ['%.7f' % i for i in evaluate_dataset(metrics, model, val_loader, device,properties)])
    if ('test' in args.split) or ('all' in args.split):
        logging.info('Testing split...')
        results.append(['test'] + ['%.7f' % i for i in evaluate_dataset(metrics, model, test_loader, device,properties)])
    header = ','.join(header)
    results = np.array(results)

    np.savetxt(os.path.join(args.modelpath, 'evaluation.csv'), results, header=header, fmt='%s', delimiter=',')


def evaluate_dataset(metrics, model, loader, device,properties):
    # TODO: Adapt for SCHNARC, copy old
    for metric in metrics:
        metric.reset()

    qm_values={}
    predicted={}
    header=[]
    for batch in loader:
        batch = {
            k: v.to(device)
            for k, v in batch.items()
        }
        result = model(batch)

        for prop in result:
            if prop in predicted:
                predicted[prop] += [result[prop].cpu().detach().numpy()]
            else:
                predicted[prop] = [result[prop].cpu().detach().numpy()]
        for prop in batch:
            if prop in qm_values:
                qm_values[prop] += [batch[prop].cpu().detach().numpy()]
            else:
                qm_values[prop] = [batch[prop].cpu().detach().numpy()]

        #for metric in metrics:
            #metric.add_batch(batch, result)
    results = [metric.aggregate() for metric in metrics]

    import numpy as np
    np.savez("prediction.npz",predicted)
    np.savez("reference.npz",qm_values)
    #logging.info('Stored model predictions in {:s} ...'.format(prediction_path))

    return results


def main(args):
    # setup
    train_args = setup_run(args)

    device = torch.device("cuda" if args.cuda else "cpu")
 
    # get dataset
    environment_provider = get_environment_provider(train_args, device=device)
    dataset = get_dataset(train_args, environment_provider=environment_provider)

    # get dataloaders
    split_path = os.path.join(args.modelpath, "split.npz")
    train_loader, val_loader, test_loader = get_loaders(
        args, dataset=dataset, split_path=split_path, logging=logging
    )
    
    # define metrics
    metrics = get_metrics(train_args)
    
   
    # train or evaluate
    if args.mode == "train":

        # get statistics
        atomref = dataset.get_atomref(args.property)
        if args.property == "charges":
            mean = {}
            stddev = {}
            mean[args.property],stddev[args.property] = None,None 
        else:
            mean, stddev = get_statistics(
                args=args,
                split_path=split_path,
                train_loader=train_loader,
                atomref=atomref,
                divide_by_atoms=get_divide_by_atoms(args),
                logging=logging,
            )

        # build model
        model = get_model(args, train_loader, mean, stddev, atomref, logging=logging)

        # build trainer
        logging.info("training...")
        trainer = get_trainer(args, model, train_loader, val_loader, metrics)

        # run training
        trainer.train(device, n_epochs=args.n_epochs)
        logging.info("...training done!")

    elif args.mode == "eval":

        # remove old evaluation files
        evaluation_fp = os.path.join(args.modelpath, "evaluation.txt")
        if os.path.exists(evaluation_fp):
            if args.overwrite:
                os.remove(evaluation_fp)
            else:
                raise ScriptError(
                    "The evaluation file does already exist at {}! Add overwrite flag"
                    " to remove.".format(evaluation_fp)
                )

        # load model
        logging.info("loading trained model...")
        model = spk.utils.load_model(
            os.path.join(args.modelpath, "best_model"), map_location=device
        )
        # run evaluation
        logging.info("evaluating...")
        
        if spk.utils.get_derivative(train_args) is None:
            with torch.no_grad():
                evaluate(args,model,train_loader,val_loader,test_loader,device,metrics=metrics,)
        else:
            evaluate(args,model,train_loader,val_loader,test_loader,device,metrics=metrics,)
        #evaluate_(args,model,train_loader,val_loader,test_loader,device,train_args)
        logging.info("... evaluation done!")

    else:

        raise ScriptError("Unknown mode: {}".format(args.mode))


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()

    if args.mode == "from_json":
        args = spk.utils.read_from_json(args.json_path)

    main(args)
