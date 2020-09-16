import io
import numpy as np
from random import shuffle
from gzip import GzipFile as gzip
from random import shuffle
from sklearn.neighbors import NearestNeighbors


mig = []
y = []

infile=np.loadtxt('input.gen')
#Remove columns with missing data.
infile=np.ma.compress_cols(np.ma.masked_invalid(infile))
num_samples=150
for i in range(0,num_samples):
	idx = np.random.choice(infile.shape[1], 130, replace=False)
	n = infile[:,idx]
	mig.append(np.array(n).T)
x = np.array(mig)



from keras.models import load_model
from sklearn.metrics import confusion_matrix

model = load_model('Trained_Model.acc.mod')
pred = model.predict(x)
print(pred)
print(np.mean(pred, axis=0))
