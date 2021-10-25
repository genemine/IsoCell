import pandas as pd
import numpy as np
from sys import argv

script,data,label,classes = argv
classes=int(classes)

exp=pd.read_csv(data)
targetid=exp.index
exp=exp.T
exp=exp.reset_index()
label=pd.read_csv(label)
exp['label']=label
exp=exp.sort_values("label")

exp=exp.groupby('label').median()
exp=exp.T
#exp.to_csv('./tmp.csv',index=False)


median=exp.reset_index(drop=True)
cellfc=pd.DataFrame()
i=1
while i<=classes :
    j=i+1
    while j<=classes :
        rowx=i
        rowy=j
        row=rowx+rowy
        cellfc[row]=(median[rowx]+1)/(median[rowy]+1)
        j=j+1
    i=i+1
#cellfc=cellfc.replace(np.inf,-1)

cellfc=cellfc.values
cellfc=np.log2(cellfc)
cellfc=np.abs(cellfc)
cellfc=pd.DataFrame(cellfc)


averfc=pd.DataFrame()
averfc['targetid']=targetid
averfc['mean_value']=cellfc.mean(axis=1)
averfc=averfc.sort_values(by="mean_value" , ascending=False)
averfc.to_csv('./fcmean.csv',sep=',',index=False)



