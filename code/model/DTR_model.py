from sklearn.model_selection import LeaveOneOut
from sklearn.tree import DecisionTreeRegressor
import numpy as np

def DTR_model(train_data,test_data,params):
    train = open(train_data,'r')
    data = train.readlines()
    x=[];y=[]
    num_sample = len(data)
    for i in range(num_sample):
        sample=data[i].split('\t')
        num_feature = len(sample)-1
        sample_feature=[]
        for n in range(0,num_feature):
            sample_feature.append(float(sample[n]))
        x.append(sample_feature)
        y.append(float(sample[-1]))
    test_sample = [test_data]
    loo = LeaveOneOut()
    y_predict=[]
    for train, test in loo.split(x):
        x_train, x_test, y_train, y_test = np.array(x)[train], np.array(x)[test], np.array(y)[train], np.array(y)[test]
        model = DecisionTreeRegressor(**params)
        model.fit(x_train, y_train)
        y_out = model.predict(test_sample)
        y_predict.append(y_out)
    result=sum(y_predict)/len(y_predict)
    return result[0]