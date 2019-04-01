import io
import numpy as np
from random import shuffle
from gzip import GzipFile as gzip
from random import shuffle
from sklearn.neighbors import NearestNeighbors

def sort_min_diff(amat):
    '''this function takes in a SNP matrix with indv on rows and returns the same matrix with indvs sorted by genetic similarity.
    this problem is NP-hard, so here we use a nearest neighbors approx.  it's not perfect, but it's fast and generally performs ok.
    assumes your input matrix is a numpy array'''
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

mig = []
y = []

infile=np.loadtxt('input.gen')
#Remove columns with missing data.
infile=np.ma.compress_cols(np.ma.masked_invalid(infile))
num_samples=150
for i in range(0,num_samples):
	idx = np.random.choice(infile.shape[1], 130, replace=False)
	n = infile[:,idx]
	#n = sort_min_diff(n[:,:])
	mig.append(np.array(n).T)
#y = [0 for i in xrange(100)]
#y.extend([1 for i in xrange(100)])
#y.extend([2 for i in xrange(100)])
#y.extend([3 for i in xrange(100)])
x = np.array(mig)

#print len(x), len(y)


from keras.models import load_model
from sklearn.metrics import confusion_matrix

model = load_model('big.data.89.2.acc.mod')
pred = model.predict(x)
print(pred)
print(np.mean(pred, axis=0))
