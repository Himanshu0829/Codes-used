#!/usr/bin/env python
# coding: utf-8

# In[46]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns


data1=pd.read_csv('plot_pop_rmsd.csv', header=None)
data2=pd.read_csv('plot_co_ma.csv', header=None)
data3=pd.read_csv('plot_int_ma.csv', header=None)
data4=pd.read_csv('plot_ma_element.csv', header=None)
data5=pd.read_csv('plot_pop.csv', header=None)
data6=pd.read_csv('plot_pop_vs_time.csv', header=None)
#data7=pd.read_csv('analysis/plot_co_rmsd.csv', header=None)



data1.head()
#data2.head()
#data3.head()
#data4.head()
#data5.head()
#data6.head()



l=len(data6[0])

x=np.arange(l)
print(l)
x

plt.figure()


#########RMSD POP##########
row=2
col=3
fig,ax=plt.subplots(row,col)



ax=plt.subplot(2,3,1)
plt.plot(data1[0], data1[1])
plt.title('RMSD_pop',fontsize=10)
plt.xlabel('time')
plt.ylabel('RMSD')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)




#ax=plt.subplot(3,3,2)
#plt.plot(data7[0], data7[1])
#plt.title('RMSD_co',fontsize=10)
#plt.xlabel('time')
#plt.ylabel('RMSD')
#ax.xaxis.set_tick_params(labelsize=6)
#ax.yaxis.set_tick_params(labelsize=6)




ax=plt.subplot(2, 3, 2)
plt.scatter(data3[0], data3[1],s=7,c="indigo")
plt.title('Interaction Matrix',fontsize=10)
plt.xlabel('Int_Expected')
plt.ylabel('Int_Calculated')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)






ax=plt.subplot(2, 3, 3)
plt.scatter(data2[0], data2[1],s=7,c="indigo")
plt.title('Co-occurrence',fontsize=10)
plt.xlabel('Co_expected')
plt.ylabel('Co-Calculated')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)


ax=plt.subplot(2, 3, 4)
plt.plot(data4[0], data4[1])

plt.plot(data4[0], data4[2])
plt.plot(data4[0], data4[3])
plt.title('Matrix Element',fontsize=10)
plt.xlabel('time')
plt.ylabel('Ma_element')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)



ax=plt.subplot(2, 3, 5)
plt.scatter(data5[0], data5[1])
#plt.plot(data5[0],data5[2])

plt.title('Abundance',fontsize=10)
plt.xlabel('Expected')
plt.ylabel('Calculated')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)


ax=plt.subplot(2, 3, 6)
plt.plot(x,data6[0])
plt.plot(x,data6[1])
plt.plot(x,data6[2])
plt.plot(x,data6[3])

plt.title('Population',fontsize=10)
plt.xlabel('time')
plt.ylabel('Population')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)


plt.suptitle("Abundance RMSD(Asymetric_co)",fontsize=15)
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9, top=0.80, wspace=0.7, hspace=0.95)
fig.savefig('test.png',dpi=200)

plt.show()



