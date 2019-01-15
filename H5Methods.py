#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 10:43:59 2018

@author: mattmiller
"""
import numpy as np
import glob
import h5py
import random
import os
import math

#DEPRECIATED: HDF hierarchy no longer like this
def create_hdf5_file(file_name, kmer_size, contig_size, organism):
    #neg_classes = ['fungi', 'arch', 'protozoa', 'virus']
    neg_classes = [organism]
    pos_contigs = 0
    neg_contigs = 0
    with h5py.File(file_name, "w") as f:
        num_features = get_num_features(kmer_size)
        pos_group = f.create_group("Girus")
        neg_group = f.create_group("Not Girus")
        pos_glob = "/rsgrps/bhurwitz/mattmiller899/girus_work/COMPILED_FEATURES_files/girus/%d_contigs/%s/*" % (contig_size, str(kmer_size))
        for input_file in glob.glob(pos_glob):
            if os.path.isdir(input_file):
                continue
            org_name = os.path.basename(input_file).split(".")[0]
            org_group = pos_group.create_group(org_name)
            with open(input_file, "r") as infile:
                line_count = 0
                for l in infile:
                    larr = np.fromstring(l, sep = ",")
                    contig_dset = org_group.create_dataset(name = "c%d" % line_count,
                                                           data = larr,
                                                           shape=(1, num_features),
                                                           dtype=np.float64,
                                                           compression='gzip')
                    line_count += 1
                    pos_contigs += 1

        for neg_class in neg_classes:
            neg_glob = "/rsgrps/bhurwitz/mattmiller899/girus_work/COMPILED_FEATURES_files/%s/%d_contigs/%s/*" % (neg_class, contig_size, str(kmer_size))
            for input_file in glob.glob(neg_glob):
                if os.path.isdir(input_file):
                    continue
                org_name = os.path.basename(input_file).split(".")[0]
                org_group = neg_group.create_group(org_name)
                with open(input_file, "r") as infile:
                    line_count = 0
                    for l in infile:
                        larr = np.fromstring(l, sep = ",")
                        contig_dset = org_group.create_dataset(name = "c%d" % line_count,
                                                               data = larr,
                                                               shape=(1, num_features),
                                                               dtype=np.float64,
                                                               compression='gzip')
                        line_count += 1
                        neg_contigs += 1
    
    return pos_contigs, neg_contigs
    

def get_num_features(kmer):
    if kmer == "no_kmer":
        return int(66)
    else:
        return int(66 + ((4 ** int(kmer)) / 2))

def get_num_features_kmer_only(kmer):
    return int((4 ** int(kmer)) / 2)
    
def printname(name, item):
    return name


def h5py_dataset_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)


def h5py_group_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Group):
            yield (path, item)
            

def create_training_validation_sets(file_name, train_pos_contigs, train_neg_contigs):
    partition = {}
    train = []
    validation = []
    labels = {}
    total_groups = []
    with h5py.File(file_name, "r") as f:
        #Do 2 main groups: Girus and Not Girus
        g_group = f.get("Girus")
        for (top_path, top_group) in h5py_group_iterator(f):
            if top_group.name == "Girus":
                ID = 1
            else:
                ID = 0
            for (path, group) in h5py_group_iterator(top_group, prefix = top_group.name):
                total_groups.append((path, group))
                """
                if random.random() <= fraction_train:
                    print("%s %s" % (path, group))
                    for (dpath, dset) in h5py_dataset_iterator(group, prefix = path):
                        if train_count >= num_train_samples:
                            break
                        train.append(dpath)
                        labels[dpath] = ID
                    
                else:
                    for (dpath, dset) in h5py_dataset_iterator(group, prefix = group.name):
                        validation.append(dpath)
                        labels[dpath] = ID
                """
        #Shuffle orgs, keep track of pos and neg contig counts, compare to passed pos/neg contig counts needed for training
        random.shuffle(total_groups)
        neg_pos_counts = [0,0]
        neg_pos_num_contigs = [train_neg_contigs, train_pos_contigs]
        for (path, group) in total_groups:
            if "Not Girus" in path:
                ID = 0
            else:
                ID = 1
            names = group.keys()
            #If every contig in this org can fit in the train set, put them all in
            if neg_pos_counts[ID] + len(names) <= neg_pos_num_contigs[ID]:
                for name in names:
                    tmp_path = "{}/{}".format(path, name)
                    train.append(tmp_path)
                    neg_pos_counts[ID] += 1
                    labels[tmp_path] = ID
            #If not but theres still room, only put some of that orgs contigs in, throw away rest
            elif neg_pos_counts[ID] < neg_pos_num_contigs[ID]:
                cutoff_count = neg_pos_num_contigs[ID] - neg_pos_counts[ID]
                for i, name in enumerate(names):
                    if i == cutoff_count:
                        break
                    tmp_path = "{}/{}".format(path, name)
                    train.append(tmp_path)
                    neg_pos_counts[ID] += 1
                    labels[tmp_path] = ID
            #No room, put into validation
            else:
                for name in names:
                    tmp_path = "{}/{}".format(path, name)
                    validation.append(tmp_path)
                    labels[tmp_path] = ID
        
        print("Train len: %d" % len(train))
        print("Validation len: %d" % len(validation))
        partition['train'] = train
        partition['validation'] = validation
        return partition, labels
    

def old_walk_hdf(hdf_file):
    partition = []
    labels = {}
    total_groups = []
    pos_neg_contigs = [0,0]
    top_levels = ["/Not Girus", "/Girus"]
    with h5py.File(hdf_file, "r") as f:
        #Do 2 main groups: Girus and Not Girus
        for i, top in enumerate(top_levels):
            top_group = f.get(top)
            for (path, group) in h5py_group_iterator(top_group):
                names = group.keys()
                for name in names:
                    tmp_path = "{}{}/{}".format(top, path, name)
                    partition.append(tmp_path)
                    labels[tmp_path] = i
                    pos_neg_contigs[i] += 1
    #return 1 then 0 for pos/neg
    return partition, labels, pos_neg_contigs[1], pos_neg_contigs[0]


def walk_hdf(hdf_file):
    partition = []
    labels = {}
    total_groups = []
    pos_neg_contigs = [0,0]
    top_levels = ["/Not Girus", "/Girus"]
    with h5py.File(hdf_file, "r") as f:
        #Do 2 main groups: Girus and Not Girus
        for i, top in enumerate(top_levels):
            top_dataset = f.get(top)
            num_contigs = len(top_dataset)
            print('len of {} = {}'.format(top, num_contigs))
            for j in range(num_contigs):
                tmp_path = '{}/{}'.format(top, j)
                partition.append(tmp_path)
                labels[tmp_path] = i
            pos_neg_contigs[i] += num_contigs
    #return 1 then 0 for pos/neg
    return partition, labels, pos_neg_contigs[1], pos_neg_contigs[0]


def walk_hdf_multiclass(hdf_file):
    partition = []
    labels = {}
    total_groups = []
    pos_neg_contigs = [0,0]
    top_levels = ["/Girus", '/Virus', '/Fungi', '/Protozoa', '/Arch', '/Bact']
    with h5py.File(hdf_file, "r") as f:
        #Do 2 main groups: Girus and Not Girus
        for i, top in enumerate(top_levels):
            top_group = f.get(top)
            for (path, group) in h5py_group_iterator(top_group):
                names = group.keys()
                for name in names:
                    tmp_path = "{}{}/{}".format(top, path, name)
                    partition.append(tmp_path)
                    labels[tmp_path] = i
                    pos_neg_contigs[i] += 1
    return partition, labels, pos_neg_contigs

def load_entire_hdf(hdf_file):
    pos_neg_contigs = [0,0]
    top_levels = ["/Not Girus", "/Girus"]
    X_arr = []
    Y_arr = []
    with h5py.File(hdf_file, "r") as f:
        #Do 2 main groups: Girus and Not Girus
        for i, top in enumerate(top_levels):
            top_dataset = f.get(top)
            num_contigs = len(top_dataset)
            X_arr.append(np.array(top_dataset))
            Y_arr.append(np.array([i for j in range(num_contigs)]))
            pos_neg_contigs[i] += num_contigs
    X = np.append(X_arr[0], X_arr[1], axis=0)
    y = np.append(Y_arr[0], Y_arr[1], axis=0)
    return X, y, pos_neg_contigs[1], pos_neg_contigs[0]


def load_pos_neg_hdf(hdf_file):
    top_levels = ["/Not Girus", "/Girus"]
    X_arr = []
    Y_arr = []
    with h5py.File(hdf_file, "r") as f:
        #Do 2 main groups: Girus and Not Girus
        for i, top in enumerate(top_levels):
            top_dataset = f.get(top)
            num_contigs = len(top_dataset)
            X_arr.append(np.array(top_dataset))
            Y_arr.append(np.array([i for j in range(num_contigs)]))
    return X_arr[0], Y_arr[0], X_arr[1], Y_arr[1]

