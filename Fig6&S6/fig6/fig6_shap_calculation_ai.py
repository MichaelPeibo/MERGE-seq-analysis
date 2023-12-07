
#Run the below code in your notebook to check the installed version
import shap
import matplotlib.pyplot as plt
from pycaret.classification import *
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
import random
plt.rcParams['pdf.fonttype'] = 'truetype'
import numpy as np
import pandas as pd

df=pd.read_csv('./shap_results/fig6_var100_AI_valid_binary_normdata_allcells.csv')
exp1 = setup(df, target = 'binary',session_id=1,index=False,n_jobs=4)
nb = create_model('nb')
Xtest = get_config('X_test')
ytest = get_config('y_test')
X_train = get_config('X_train')
y_train = get_config('y_train')
sample_size=1500
sub_sampled_train_data = shap.sample(X_train, sample_size, random_state=0)
sub_sampled_test_data = shap.sample(Xtest, sample_size, random_state=0)
explainer = shap.KernelExplainer(nb.predict, sub_sampled_train_data)## nb,predict_proba generates barplot
shap_values = explainer.shap_values(sub_sampled_test_data)
np.save("./shap_results/fig6_ai_shap_values_1500cells.npy", shap_values)
shap.summary_plot(shap_values, sub_sampled_test_data,show=False)
plt.savefig("./shap_results/fig6_ai var shap plot nb 1500cells.pdf",dpi=1000)