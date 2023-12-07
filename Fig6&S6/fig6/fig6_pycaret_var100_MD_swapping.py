
import matplotlib.pyplot as plt
from pycaret.classification import *
from sklearn.utils import shuffle
import random
import pandas as pd
plt.rcParams['pdf.fonttype'] = 'truetype'
sample_trials=list(range(1,101))
for x in sample_trials:
        df = pd.read_csv('./swapping_results/fig6_var100_MD_valid_swapping_%d_1000cells.csv' % (x))
        AUC = []
        Accuracy = []
        Recall = []
        Prec = []
        F1 = []
        Kappa = []
        MCC = []
        
        setup_data = setup(df, target='binary', session_id=1, index=False, n_jobs=4)
        nb = create_model('nb')
        Xtest = get_config('X_test')
        ytest = get_config('y_test')
        predict_model(nb)
        
        AUC.append(pull()['AUC'][0])
        Accuracy.append(pull()['Accuracy'][0])
        Recall.append(pull()['Recall'][0])
        Prec.append(pull()['Prec.'][0])
        F1.append(pull()['F1'][0])
        Kappa.append(pull()['Kappa'][0])
        MCC.append(pull()['MCC'][0])
        
        metrics_data = pd.DataFrame({
            'AUC': AUC,
            'Accuracy': Accuracy,
            'Recall': Recall,
            'Prec.': Prec,
            'F1': F1,
            'Kappa': Kappa,
            'MCC': MCC
        })
        
        metrics_data.to_csv(f"./MD_results_swapping_1000cells/MD_var_model_results_genes_metrics_{x}.csv", index=False)
