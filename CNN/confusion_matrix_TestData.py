import numpy as np
import tensorflow as tf
from random import shuffle
from gzip import GzipFile as gzip
import io
import numpy as np
from sklearn.neighbors import NearestNeighbors

u1 = np.load("simModel1.npz")
u2 = np.load("simModel2.npz")
u3 = np.load("simModel3.npz")
x=np.concatenate((u1['simModel1'],u2['simModel2'],u3['simModel3']),axis=0)

y=[0 for i in xrange(len(u1['simModel1']))]
y.extend([1 for i in xrange(len(u2['simModel2']))])
y.extend([2 for i in xrange(len(u3['simModel3']))])


print len(x), len(y)


from keras.models import load_model
from sklearn.metrics import confusion_matrix

model = load_model('Trained_Model.acc.mod')
pred = model.predict(x)
pred_cat = [i.argmax() for i in pred]
print confusion_matrix(y, pred_cat)
print
print confusion_matrix(y, pred_cat) / float(len(y))

k = []
for idx,i in enumerate(pred):
    k.append(i/sum(i))

n = []
for i,j in zip(k,y):
    if i.argmax() == j: val = 1
    else: val=0
    prob = i[i.argmax()]
    n.append((prob, val))

d = {}
s,e = 0,.1
for i in range(10):
    if e > .999: e=1
    d[(s,e)] = [0., 0.]
    s+=.1
    e+=.1

for i,o in n:
    for j1,j2 in d:
        if j1 < i <= j2:
            d[(j1,j2)][o]+=1

for i in sorted(d):
    if sum(d[i]): print i,'\t', d[i][1] / sum(d[i]), '\t', sum(d[i])
