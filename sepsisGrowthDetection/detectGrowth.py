'''
Created on Feb 25, 2015

@author: michaelwalton
'''
from os import path
from os import listdir
from os import chdir
import numpy as np
from datetime import datetime
from sklearn.cross_validation import train_test_split
from itertools import izip
import cPickle
import time

#input: path to metadata file
#2dim: array of strings containing raw metadata
def loadMetadata(mdpath, ncol=None):
    if(ncol==None):
        md=np.genfromtxt(mdpath, delimiter=",", dtype=None)
    else:
        md=np.genfromtxt(mdpath, delimiter=",", dtype=None, usecols=np.arange(0,ncol))
        
    #trim off the header
    return md[1:,:]

#input: path to a single spotfile
#output: a tuple containing a list of timestamps and a 2dim array of shape (n_samples, n_indicatorChannels)
def loadSpots(spath, label):
    date_conv = lambda x: datetime.strptime(x, '%Y-%m-%d_%H-%M-%S')
    
    timestamps=np.genfromtxt(spath, delimiter="\t", dtype="datetime64", usecols=[0], converters={0: date_conv})
    spotValues=np.genfromtxt(spath, delimiter="\t", dtype="float32", usecols=np.arange(1,241))
    
    #generate a copy of the specimen label for each sample
    labels=np.repeat(label, np.shape(spotValues)[0])
    
    return (timestamps, spotValues, labels)

#input a 2d array containing metadata
#output a list of paths to spotfiles
def getSpotFiles(metadata):
    #the last column in the metadata file contains the paths we want
    paths=metadata[:,-1]
    
    # replace windows style path delimiters w/ unix style
    # strip the first one so we are looking at a sub dir of working directory
    frameFolders=map(lambda s: s.replace("\\","/")[1:], paths)
    
    #list the directory, sort in lexicographic order, latest revision is the last element of this list
    revFolders=map(lambda s: sorted(listdir(s))[-1], frameFolders)
    
    #build the complete path
    fullPaths=map(lambda s1,s2: path.join(s1,s2,"spots.txt"), frameFolders, revFolders)
    
    return fullPaths
        
def binarize(a, trueCriteria):
    return (a == trueCriteria).astype(int)

def loadDataset(uids, spotPaths, labels):
    #declare empty arrays to hold all the data
    times = np.zeros((1))
    x_array = np.zeros((1,240))
    y_array = np.zeros((1))
    
    #the main body of the dataset is contained in the 1d y_array and the 2d x_array
    #uid_idx is a dict mapping UIDs to indicies of the x and y arrays to separate by
    #specimen
    uid_idx = {}
    acc=0
    
    #iterate over the uids, spotfile pointers, labels and build the dataset
    for uid, xP, l in izip(uids, spotPaths, labels):
        t, x, y = loadSpots(xP, l)
        if (x.size > 240):
            uid_idx[uid] = acc
            times=np.concatenate((times, t))
            x_array=np.concatenate((x_array, x), axis=0)
            y_array=np.concatenate((y_array, y))
            acc += np.shape(x_array)[0]
            
    #trim off the first row of each array
    times=times[1:]
    y_array=y_array[1:]
    x_array=x_array[1:]
    
    return (times, x_array, y_array, uid_idx)
    
#input: the metadata file
def buildDataset(metadata):
    #Because the sample order matters in the test set, we want to separate this from the complete
    #dataset before we start aggregating the training set (which will be mixed across specimens)
    uids=metadata[:,0]
    specimenLabels=metadata[:,4]
    labels = binarize(specimenLabels, "Control")
    spotPaths=getSpotFiles(metadata)
    
    #this simplistic, I may need to find a split that considers the class imbalance
    id_train, id_test, x_trainPaths, x_testPaths, y_trainBase, y_testBase = train_test_split(
                                                                                             uids, 
                                                                                             spotPaths, 
                                                                                             labels, 
                                                                                             test_size=0.25,
                                                                                             random_state=123)
    
    #'dereference' the dataset file pointers and aggregate everything into contiguous blocks
    t_train, x_train, y_train, uid_train = loadDataset(id_train, x_trainPaths, y_trainBase)
    t_test, x_test, y_test, uid_test = loadDataset(id_test, x_testPaths, y_testBase)
    
    dataset = {'t_train': t_train,
               'x_train': x_train,
               'y_train':y_train,
               'uid_train':uid_train,
               't_test':t_test, 
               'x_test':x_test,
               'y_test':y_test,
               'uid_test':uid_test
               }
    
    return dataset
    
if __name__ == '__main__':
    dataPath="/Users/michaelwalton/workspace/2014-04--Sepsis-TF"
    chdir(dataPath)
    
    metaDataFile="2015-02-16_14-16-17_2014-04--Sepsis-TF_YC.csv"
    metadata=loadMetadata(metaDataFile, 7)

    #if the serialized dataset has already been generated, load that
    #else go get a cup of coffee 'cause its gonna be a while
    if (path.isfile('dataset.pkl')):
        print "Loading pickled dataset"
        start=time.clock()
        
        with open('dataset.pkl', 'rb') as handle:
            dataset = cPickle.load(handle)
        
        print("Load from pickle %s seconds" % str(time.clock() - start))
    else:
        print "Loading dataset from metadata file"
        
        start=time.clock()
        
        dataset = buildDataset(metadata)
        
        print("Load from file took %s seconds" % str(time.clock() - start))
        start=time.clock()
    
        with open('dataset.pkl','wb') as handle:
            cPickle.dump(dataset, handle)
            
        print("Pickling took %s seconds" % str(time.clock() - start))
        
    
    
    
    