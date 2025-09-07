import logging
logger = logging.getLogger(__name__)

def run_ml(**kwargs):
    """
    Auto-wrapped from 5_Machine_Learning.py.
    kwargs can pass configuration and file paths.
    """
    #read data
    input_data = "ready_for_ML_processed.csv"
    GSEA = "GOZGIT_ESR1_TARGETS_UP_GSEA.csv"
    scores = "GOZGIT_ESR1_TARGETS_UP_scores-afterML.csv"
    output = "GOZGIT_ESR1_TARGETS_UP-afterML.csv"

    #function to sep class1 from class 0 for training set
    def categorize(df, threshold):
        row = df['GSEA_score']
        df['class'] = np.where( row > threshold, 1, 0)

    #import all packages
    import warnings
    import pandas as pd
    import numpy as np
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    from sklearn import metrics
    from sklearn.utils.class_weight import compute_class_weight, compute_sample_weight
    from sklearn.model_selection import cross_validate
    from sklearn.model_selection import cross_val_predict
    from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
    from sklearn.metrics import plot_roc_curve

    #turn off warnings
    warnings.filterwarnings('ignore')

    #read in structural data
    data = pd.read_csv(input_data)

    #read in GSEA data
    scores = pd.read_csv(GSEA)

    #add GSEA info to main df
    data["GSEA_score"] = "Unknown"

    for i in range(len(data["cell_line"])):
        if data["cell_line"][i] in scores.index:
            data['GSEA_score'][i] = scores.loc[data["cell_line"][i]][0]

    #assign preliminary class1 or class0
    categorize(data, up_thresh)

    #balance data as needed by oversampling
    class1 = new_data[new_data['class']==1]
    class0 = new_data[new_data['class']==0]

    print(len(class1),len(class0))

    #perform machine learning - this is on a per protein basis 
    X = new_data_b[['x', 'y', 'z','cluster','dist','num_members']]
    y = new_data_b['class']

    rf_r = RandomForestClassifier(n_estimators=300, random_state=42, class_weight= 'balanced')

    scoring_choice = {'accuracy' : make_scorer(accuracy_score), 
                              'precision' : make_scorer(precision_score),
                              'recall' : make_scorer(recall_score), 
                              'f1_score' : make_scorer(f1_score),
                              'mcc:' : make_scorer(matthews_corrcoef)}

    rf_r_scores = cross_validate(estimator = rf_r, X=X, y=y, cv=5,
                                          scoring = scoring_choice,
                                          return_train_score=True)
    rf_r_results = cross_val_predict(rf_r, X, y, cv=5)

    #merge machine learning output with previous df
    df_results = X.copy()
    df_results['class'] = pd.Series(rf_r_results)
    new_data_final = pd.merge(df_results,new_data_b,on=['x','y','z'])

    #save ML output and scores
    rf_r_scores.to_csv(scores,index=False)
    new_data_final.to_csv(output, index=False)
    logger.info('Completed run_ml.')
    return True
