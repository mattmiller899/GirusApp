#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 10:40:24 2018

@author: mattmiller
"""
import keras
import numpy as np
import h5py
from H5Methods import *
import time
from pympler.tracker import SummaryTracker


class DataGenerator(keras.utils.Sequence):
    def __init__(self, list_IDs, labels, kmer_size, file_name, minibatch = 64, shuffle = True, tracker = None):
        self.list_IDs = list_IDs
        self.labels = labels
        self.file_name = file_name
        self.minibatch = minibatch
        self.shuffle = shuffle
        self.num_features = get_num_features(kmer_size)
        #self.tracker = tracker
        self.on_epoch_end()
    
    
    def __getitem__(self, index):
        if (index + 1) * self.minibatch > len(self.list_IDs):
            indexes = self.indexes[index * self.minibatch : len(self.list_IDs)]
            batch_size = len(self.list_IDs) - (index * self.minibatch)
        else:
            indexes = self.indexes[index * self.minibatch : (index + 1) * self.minibatch]
            batch_size = self.minibatch
        list_IDs_temp = [self.list_IDs[k] for k in indexes]
        X, y = self.__data_generation(list_IDs_temp, batch_size)
        return X, y
   
    
    def get_len(self):
        print(int(np.ceil(len(self.list_IDs) / self.minibatch)))  
        #print(int(np.floor(len(self.list_IDs) / self.minibatch)))
    
    def __len__(self):
        return int(np.ceil(len(self.list_IDs) / self.minibatch))
        #return int(np.floor(len(self.list_IDs) / self.minibatch))
        
    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffle == True:
            np.random.shuffle(self.indexes)
            
    
    def __data_generation(self, list_IDs_temp, batch_size):
        x = np.empty((batch_size, self.num_features))
        y = np.empty((batch_size))
        with h5py.File(self.file_name, "r") as f:
            pos_dataset = f.get('/Girus')
            neg_dataset = f.get('/Not Girus')
            for i, ID in enumerate(list_IDs_temp):
                ID_split = ID.split('/')
                dataset_name = ID_split[1]
                ID_index = int(ID_split[2])
                if 'Not Girus' in dataset_name:
                    x[i,] = neg_dataset[ID_index]
                else:
                    x[i,] = pos_dataset[ID_index]
                y[i] = self.labels[ID]
            #print("data generation summary")
            #self.tracker.print_diff()
            return x, y
        
        
