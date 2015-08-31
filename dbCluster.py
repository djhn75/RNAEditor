'''
Created on 13.02.2015

@author: david
'''
import numpy as np
from random import shuffle
from sklearn.metrics.pairwise import pairwise_distances 


class dbCluster():
    def __init__(self):
        pass


    
    def dbscan(self,X, eps=10, min_samples=5):
        """Perform DBSCAN clustering from vector array.
    
        Parameters
        ----------
        X: array [int1,int1]
            Array of Samples. 
            In this case it should be the positions of the variations in the genome
    
        eps: float, optional
            The maximum distance between two samples for them to be considered
            as in the same neighborhood.
    
        min_samples: int, optional
            The number of samples in a neighborhood for a point to be considered
            as a core point.
    
    
        Returns
        -------
        core_samples: array [n_core_samples]
            Indices of core samples.
    
        labels : array [n_samples]
            Cluster labels for each point.  Noisy samples are given the label -1.
    
        """
        if not eps > 0.0:
            raise ValueError("eps must be positive.")
    
        X = np.asarray(X)
        n = X.shape[0] #get number of elements (not sure) 
    
        index_order=range(n)
        shuffle(index_order)
        
        
        distance_matrix = True
        D = self.calculate1dDistanceMatrix(X,eps=eps)
    
        # Calculate neighborhood for all samples. This leaves the original point
        # in, which needs to be considered later (i.e. point i is the
        # neighborhood of point i. While True, its useless information)
        neighborhoods = []
    
        #neighborhoods = [np.where(x <= eps)[0] for x in D]
        neighborhoods=D
        # Initially, all samples are noise.
        labels = -np.ones(n, dtype=np.int)
    
        # A list of all core samples found.
        core_samples = []
    
        # label_num is the label given to the new cluster
        label_num = 0
    
        # Look at all samples and determine if they are core.
        # If they are then build a new cluster from them.
        for index in index_order:
            # Already classified
            if labels[index] != -1:
                continue
    
            # get neighbors from neighborhoods or ballTree
            index_neighborhood = []
    
            index_neighborhood = neighborhoods[index]
            
    
            # Too few samples to be core
            if len(index_neighborhood) < min_samples:
                continue
    
            core_samples.append(index)
            labels[index] = label_num
            # candidates for new core samples in the cluster.
            candidates = [index]
    
            while len(candidates) > 0:
                new_candidates = []
                # A candidate is a core point in the current cluster that has
                # not yet been used to expand the current cluster.
                for c in candidates:
                    c_neighborhood = []
                    
                    c_neighborhood = neighborhoods[c]
                    
                    noise = np.where(labels[c_neighborhood] == -1)[0] #indexes of candidate neigbours which do not belong to a cluster yet
                    noise = c_neighborhood[noise]
                    labels[noise] = label_num
                    for neighbor in noise:
                        n_neighborhood = []
                        
                        n_neighborhood = neighborhoods[neighbor]
                        
                        # check if its a core point as well
                        if len(n_neighborhood) >= min_samples:
                            # is new core point
                            new_candidates.append(neighbor)
                            core_samples.append(neighbor)
                # Update candidates for next round of cluster expansion.
                candidates = new_candidates
            # Current cluster finished.
            # Next core point found will start a new cluster.
            label_num += 1
        #return core_samples, labels
        self.coreSamples = core_samples
        self.labels = labels
    
    def calculate1dDistanceMatrix(self,lst,eps=20):
        '''
        creates a distance matrix for the given vector
        
        :param lst: vector of samples
        
        :return: np.array(diffMatrix)
        '''
        
        if not isinstance(lst, (list, tuple, np.ndarray)):
            raise TypeError("Paramer has to be eithe a List or a Tuple found %s" % type(lst))
        if not all(isinstance(item, (int,float)) for item in lst):
            raise TypeError("List should only contain numbers")
        lst = np.asarray(lst)
        diffMatrix=[]
        i = 0
        for l1 in lst:
            diffList=[]
            
            diffList= abs(lst-l1)
            diffList = np.where(diffList<=eps)[0]
            diffMatrix.append(diffList)
            
            """i+=1
            
            if i % 1000 == 0:
                print "%d distances lists calculated" % i"""
                
        return np.asarray(diffMatrix)
    
    """def getSilhuetteScore(self,X,labelSet):
        return np.mean(silhouette_samples(X, labelSet))
    """