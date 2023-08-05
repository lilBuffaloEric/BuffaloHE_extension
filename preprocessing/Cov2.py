import numpy as np
import pandas as pd
from scipy.special import expit
import matplotlib.pyplot as plt
import numpy.random as rnd
from sklearn.metrics import r2_score
from PIL import Image
import read_data


rows = 200
cols = 256
path = "EV\\EigenVectors"
DominantVecsNum = 4
Iternum = 14

# X_ori = read_data.parkinson_data()
# X_ori = read_data.air_quality_data()
# X_ori = read_data.load_yale()
# X_ori = read_data.wine_data("red")
X_ori = read_data.load_mnist(no_of_imgs=rows)  # original data matrix rows x cols
# X_ori = read_data.load_fashion_mnist(no_of_imgs=rows)


# TO GENERATE TAILORED DATASET IN .CSV FILE FOR HOMOMORPHIC PCA, PLEASE UNCOMMENT THESE TWO LINES:
#X_ori_df = pd.DataFrame(X_ori) #  columns=['col'+str(i) for i in range(1, cols+1)]  #  index=['row'+str(i) for i in range(1, rows+1)]
#X_ori_df.insert(0,'idx',['row'+str(i) for i in range(1, rows+1)])

# TO GENERATE TAILORED DATASET IN .CSV FILE FOR HOMOMORPHIC PCA, UNCOMMENT ONLY THE LINE CORRESPONDING TO THE TAILORED DATASET OBTAINED BY the read_data funtion invoked above:
# X_ori_df.to_csv(path_or_buf="mnist.csv",index=False)
# X_ori_df.to_csv(path_or_buf="fashion_mnist.csv",index=False)
#X_ori_df.to_csv(path_or_buf="yale.csv",index=False)
# X_ori_df.to_csv(path_or_buf="pks.csv",index=False)
# X_ori_df.to_csv(path_or_buf="wineRed.csv",index=False)


X_oriT = np.transpose(X_ori)
# input_matrix = (input_matrix / 255.0) * 0.99 + 0.01 # do the Scaling.
mu = np.mean(X_ori, axis=0)  # mean vector.
aggregate = np.sum(X_ori,axis=0)
print("Mean ", mu.shape)
X_cent = X_ori - mu  # decentralized data matrix
print("Data after subtracting mean ", X_cent.shape, "\n")
cov = np.cov(X_cent.T)  # straightly compute the covariance matrix by numpy function
cov = np.round(cov, 5)
print("Covariance matrix ", cov.shape, "\n")

# do the Covariance Matrix ourselves:
# compute the Mean^T*Mean

meanmatrix = mu.reshape(cols, 1) @ mu.reshape(1, cols)
# compute X^TX
XTX_NonScale = (np.transpose(X_ori) @ X_ori)
XTX_Scale = XTX_NonScale * (1 / rows)
# compute Cov = 1/N * X^TX - Mean^TMean
Cov = XTX_Scale - meanmatrix
# compute Cov = 1/N * (X-meanVec)^T(X-meanVec)
Cov2 = np.transpose(X_cent) @ X_cent * (1 / rows)

# normal way computing eigenvectors and eigenvalue.
eig_val, eig_vec = np.linalg.eigh(Cov)  # np.linalg.eig(Cov)
idx = eig_val.argsort()[::-1]
eigen_values = eig_val[idx]
eigen_vectors = eig_vec[:,idx]

# Sort eigen values and corresponding eigen vectors in descending order
'''
indices = np.arange(0, len(eig_val), 1)
indices = ([x for _, x in sorted(zip(eig_val, indices))])[::-1]
eig_val = eig_val[indices]
eig_vec = eig_vec[:, indices]
print("Sorted Eigen vectors ", eig_vec)
print("Sorted Eigen values ", eig_val, "\n")

for k in range (len(eig_vec)) :
    standard_eig_vec = eig_vec[k]
    standard_eig_val = np.inner(Cov @ standard_eig_vec, standard_eig_vec) / np.inner(standard_eig_vec, standard_eig_vec)
    print("Standard ", k,"th dominant value: ",standard_eig_val)
'''

largest_eigenvalues = eigen_values[:DominantVecsNum]
largest_eigenvectors = eigen_vectors[:, :DominantVecsNum]

print("first", DominantVecsNum, "Dominant eigenValues:", largest_eigenvalues)
print("Corresponding eigenVectors:")
print(largest_eigenvectors)

# Power method computing approximate eigenvectors and eigenvalues.
# scale =0.0005
dominant_eigmtx = None
dominant_eigvec = np.random.random(cols)
dominant_eigval = None

for k in range(DominantVecsNum):
    for i in range(Iternum):
        b = np.inner(dominant_eigvec,dominant_eigvec)
        dominant_eigvec = (Cov @ dominant_eigvec)
        a = np.inner(dominant_eigvec, dominant_eigvec)
        inva = np.sqrt(1/a)
        dominant_eigvec = dominant_eigvec * inva  # normalisation.
        j = 1
        # dominant_eigvec *= scale
        # scale *= 0.1
    dominant_eigval = np.inner(Cov @ dominant_eigvec, dominant_eigvec) / np.inner(dominant_eigvec,dominant_eigvec)
    # dominant_eigvec = dominant_eigvec / np.sqrt(np.inner(dominant_eigvec, dominant_eigvec))
    print("the ", k,"th dominant value: ",dominant_eigval)
    dominant_covfactor = dominant_eigval * (dominant_eigvec.reshape(cols,1)  @ dominant_eigvec.reshape(1,cols))
    Cov = Cov - dominant_covfactor
    if k == 0 :
        dominant_eigmtx = np.copy(dominant_eigvec)
    else :
        dominant_eigmtx = np.column_stack((dominant_eigmtx, np.copy(dominant_eigvec))  )
    dominant_eigvec = np.random.random(cols)



columns  = largest_eigenvectors.shape[1]
for i in range(columns):
    norm = np.linalg.norm(largest_eigenvectors[:,i])
    largest_eigenvectors[:,i] = largest_eigenvectors[:,i]/ norm
    if np.dot(largest_eigenvectors[:,i],dominant_eigmtx[:,i])<0:
        dominant_eigmtx[:,i] = -dominant_eigmtx[:,i]

scoreV = r2_score(largest_eigenvectors,dominant_eigmtx)
X = X_ori - mu
X_red = X.dot(dominant_eigmtx)
X_new = X_red.dot(dominant_eigmtx.T) + mu
scoreX = r2_score(X_ori,X_new)



# TO TEST THE RESULT(EIGENVECTORS) OF THE HOMOMORPHIC PCA, PLEASE SET THE VARIABLE "path" AS THE DIRECTORY OF THE RESULT, AND UNCOMMENT THE FOLLOWING LINES.
'''
eigmtx_enc_df = pd.read_csv(path,header=None)
dominant_eigmtx_enc = np.array(eigmtx_enc_df).T
for i in range(columns):
    if np.dot(largest_eigenvectors[:,i],dominant_eigmtx_enc[:,i])<0:
        dominant_eigmtx_enc[:,i] = -dominant_eigmtx_enc[:,i]

scoreVE = r2_score(largest_eigenvectors,dominant_eigmtx_enc)
X = X_ori - mu
X_red = X.dot(dominant_eigmtx_enc)
X_new = X_red.dot(dominant_eigmtx_enc.T) + mu
scoreXE = r2_score(X_ori,X_new)
'''



i = 1