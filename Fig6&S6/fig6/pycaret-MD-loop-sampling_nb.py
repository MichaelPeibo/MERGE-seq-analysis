
import matplotlib.pyplot as plt
from pycaret.classification import *
from sklearn.utils import shuffle
import random
import pandas as pd
plt.rcParams['pdf.fonttype'] = 'truetype'
sample_trials=list(range(1,101))
num_vars=[2,5,10,20,50,100,200,300,400,500,1000,2000,5000]
for x in sample_trials:
    df = pd.read_csv('./sampling_results/var5000_MD_valid_binary_normdata_1000cells_%d.csv' % (x))
    for n in num_vars:
        AUC = []
        Accuracy = []
        Recall = []
        Prec = []
        F1 = []
        Kappa = []
        MCC = []
        
        input_data = df.iloc[:, 0:n]
        input_data['binary'] = df['binary']
        setup_data = setup(input_data, target='binary', session_id=1, index=False, n_jobs=4)
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
        
        metrics_data.to_csv(f"./md_results_nb/MD_var_model_results_{n}_genes_metrics_{x}.csv", index=False)
