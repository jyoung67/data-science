import numpy as np
import logging

class KnnClassifier(object):
    @staticmethod   
    def euclidean(x,y):
        """
        Inputput
            x: a 1-dimensional numpy array
            y: a 1-dimensional numpy array
        Output
            A single number, the euclidean distance between x and y
        """      

        return(np.sqrt(np.sum((x - y)**2)))
    
    @staticmethod
    def many_distances(x,Y):

        """
        Input
            x: a single point as a numpy array
            Y: a 2-dimensional numpy array
        Output
            a numpy array of euclidean distances
        """
        result = np.zeros(Y.shape[0])
        for idx in range(0, Y.shape[0]):            
            dist = KnnClassifier.euclidean(x, Y[idx])          
            result[idx] = dist  
        return(result)
    
    @staticmethod
    def closest_indices(dists, n):

        """
        Input
            dists: a numpy array of distances (1 dimensional)
            n: the number of lowest values to find
        Output
            a numpy array with the indexes in dists where the
            n lowest values are located.
        """

        arrIndexes = np.argsort(dists)   
        return(arrIndexes[0:n])
    
    @staticmethod
    def get_values_by_index(indices, values):
        """
        Input
            indices: a numpy array of indices
            values: a list of values
        Output
            a list of elements from values whose indexes
            are the values in indices
        """
        arr = np.array(values)
        return(arr[indices])
    
    @staticmethod
    def get_mode(values):
        """
        Input
            values: a lists of values
        Output
            the most common value in the list.
            If there's a tie, break it any way you want to.
        """

        counts = np.unique(values, return_counts=True)
        cnt = np.argsort(counts[1])   
        colorIndex = cnt[len(cnt) -1]
        label = counts[0][colorIndex]   
        return(label)
    
    @staticmethod
    def knn(ind, test_pts, train_pts, labels, k):
        """
        Input
            test: A 2-D numpy array of points (rows) and features (cols)
            train: A 2-D numpy array of points (rows) and features (cols)
            labels: A list of labels associated with train points
        Output
            A list of best guesses for the test labels
        """      
        output = []    
        startIdx = ind[0]
        endIdx = ind[1]
        subsetOfPoints = test_pts[startIdx:endIdx]
       
        for i in range(0,subsetOfPoints.shape[0]):                     
            dists = KnnClassifier.many_distances(subsetOfPoints[i],train_pts)              
            label_indices = KnnClassifier.closest_indices(dists, k)           
            labelSet = KnnClassifier.get_values_by_index(label_indices, labels)          
            predictedLabel = KnnClassifier.get_mode(labelSet)
            output.append(predictedLabel)      
            
        
        return output
    
    @staticmethod    
    def multiclass_accuracy(truth,predictions):
        """
        Input
        truth: a list of true labels
        predictions: a list of predicted labels
        Output
        a single value - the multiclass accuracy
        """
        cmp = np.array(truth) == np.array(predictions)
        wrong = len(cmp[cmp == False])
        right = len(cmp[cmp == True])
        accuracy = right/(wrong + right)
        return(accuracy)