#!/usr/bin/env python
import pandas as pd
import argparse
import sys
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
#import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description = 'linear model to predict sera-Ag reaction')
parser.add_argument('-a', '--alpha', type = float, default = 0.5, help = 'higher to reduce model complexity; default 0.5')
parser.add_argument('-e', '--elisa', help = 'file with elisa readings' )
parser.add_argument('-p', '--pair_diff', help = 'pair diff file')
args = parser.parse_args()

# read files
elisa = pd.read_csv(args.elisa, sep=',', header = 'infer', names = None, index_col = None)
seq_diff = pd.read_csv(args.pair_diff, sep='\t', header = 'infer', names = None, index_col = None)
#print(elisa.shape)
#print(elisa.info())
#print(seq_diff)
#print(seq_diff.info())
#sys.exit()

def change_col_name(row):
    return row['ab'].replace('anti-', '')

def add_seq_diff(row):
    if row['antigen'] == row['ab2']:
        return pd.Series([0.0] * (len(diff_wide.columns)-2), index = windows)
    else:
        ab1 = row['antigen']
        ab2 = row['ab2']
        selected = diff_wide[((diff_wide.id1 == ab1) & (diff_wide.id2 == ab2)) | ((diff_wide.id2 == ab1) & (diff_wide.id1 == ab2))]
        if len(selected.index) == 1:
            return selected.iloc[0,2:len(diff_wide.columns)]
        
# filter and join data set: only native antigens and no cocktail
elisa_native_ag = elisa[elisa["antigen"].str.contains("^[A-Z]$")] # native antigens only
elisa_native_ag = elisa_native_ag.assign(ab2 = elisa_native_ag.apply(lambda x: change_col_name(x), axis = 1)) # add a column
elisa_clean = elisa_native_ag[(elisa_native_ag.ab2 != 'CKT') & ((elisa_native_ag.mouse == '1') | (elisa_native_ag.mouse == '2'))] # use blot or Ivano data for training are bad
diff_wide = seq_diff.pivot_table(index = ["id1", "id2"], # pivot to wide
                                 columns = "window",
                                 values = "diff_hydro").reset_index() # seq diff weighted by hydrophobicity

windows = diff_wide.columns[2:] # necessary for filtered windows
print(windows)
print(elisa_clean)
print(diff_wide)
diff_extracted = elisa_clean.apply(lambda x: add_seq_diff(x), axis = 1)
print(diff_extracted)
joined = elisa_clean.join(diff_extracted)
print(joined)
#sys.exit()
# Seperate into train and test
train = joined[joined['ab'].str.contains("^anti-[A-Z]$")] # 16 nat ab against 17 nat ag
test = joined[joined['ab'].str.contains("^anti-CT")] # 2 syn ab vs 17 nat ag (no cocktail)
train_pos = train[train['OD'] > 0]
test_pos = test[test['OD'] > 0]
y_train = -np.log(train_pos['OD']) # logirthm is better
y_test = -np.log(test_pos['OD'])
#y_train = train['OD']
#y_test = test['OD']
X_train = train_pos.iloc[:,8:]
X_test = test_pos.iloc[:,8:]
test_ag = test_pos['antigen']
test_ab = test_pos['antibody']
print(y_train)
print(y_test)
print(X_train)
print(X_test)
print(test_ag)
print(test_ab)
#joined.to_csv("joined.csv")
# Create linear regression object
alpha = args.alpha
#reg_ols = linear_model.LinearRegression() # really bad
reg_ridge = linear_model.Ridge(alpha = alpha) # okay
#reg_lasso = linear_model.Lasso(alpha = alpha) # bad, uniform
#reg_enet = linear_model.ElasticNet(alpha = alpha, l1_ratio=0.7) # bad, uniform
#reg_lars = linear_model.LassoLars(alpha = alpha) # bad, uniform
#reg_rob = linear_model.RANSACRegressor() # bad (no penalization)

# Train the model using the training sets
#reg_ols.fit(X_train, y_train)
reg_ridge.fit(X_train, y_train)
#reg_lasso.fit(X_train, y_train)
#reg_enet.fit(X_train, y_train)
#reg_lars.fit(X_train, y_train)
#reg_rob.fit(X_train, y_train) 

# Make predictions using the testing set
#y_ols = reg_ols.predict(X_test)
y_ridge = reg_ridge.predict(X_test)
#y_lasso = reg_lasso.predict(X_test)
#y_enet = reg_enet.predict(X_test)
#y_lars = reg_lars.predict(X_test)
#y_rob = reg_rob.predict(X_test)
#y_df = pd.DataFrame({"test": y_test, "ols": y_ols, "ridge": y_ridge, 'lasso': y_lasso, "enet": y_enet, 'lars': y_lars, 'ag': test_ag, 'ab': test_ab, 'robust': y_rob})
y_df = pd.DataFrame({"test": y_test, "ridge": y_ridge, 'ag': test_ag, 'ab': test_ab, 'od.obs': np.exp(-y_test), 'pd.pred': np.exp(-y_ridge)})
print('Mean squared error:\t%.2f' % mean_squared_error(y_test, y_ridge))
print('Coefficient of determination:\t%.2f' % r2_score(y_test, y_ridge))
# The coefficients
#print('Coefficients: \n', regr.coef_)
# The mean squared error
#print('Mean squared error:\nOLS\t%.2f\nRidge\t%.2f\nLasso\t%.2f\nEnet\t%.2f\nLars\t%.2f'
#            % (#mean_squared_error(y_test, y_ols),
#               mean_squared_error(y_test, y_ridge),
               #mean_squared_error(y_test, y_lasso),
               #mean_squared_error(y_test, y_enet),
               #mean_squared_error(y_test, y_lars)
#            )
#            )
# The coefficient of determination: 1 is perfect prediction
#print('Coefficient of determination:\nOLS\t%.2f\nRidge\t%.2f\nLasso\t%.2f\nEnet\t%.2f\nLars\t%.2f'
#            % (#r2_score(y_test, y_ols),
#               r2_score(y_test, y_ridge),
               #r2_score(y_test, y_lasso),
               #r2_score(y_test, y_enet),
               #r2_score(y_test, y_lars)
#            )
#)

pd.Series(reg_ridge.coef_).to_csv("paramter.csv")
y_df.to_excel("pred.xlsx")
'''
# Plot outputs
plt.scatter(diabetes_X_test, diabetes_y_test,  color='black')
plt.plot(diabetes_X_test, diabetes_y_pred, color='blue', linewidth=3)

plt.xticks(())
plt.yticks(())

plt.show()
'''
sys.exit()
