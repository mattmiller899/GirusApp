import logging
from keras.layers import *
import numpy as np
from sklearn.metrics import roc_auc_score
from DataGenerator import DataGenerator
from H5Methods import *
import argparse
from sklearn.metrics import average_precision_score
import time
from tensorflow.python.client import device_lib
import tensorflow as tf
import keras
import os
import resource
from keras.callbacks import Callback
import sys
import gc
from keras import backend as K
from keras.models import Sequential
from keras.layers import Dense, BatchNormalization, Dropout, Activation

def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--hdf', required=True, dest="train_file",
                            help='path to HDF file containing feature data')
    arg_parser.add_argument('-k', '--kmer', required=True, type=int,
                           help='size of kmer (can be no_kmer)')
    arg_parser.add_argument('-w', '--workers', required=True, type=int,
                            help='number of cores')
    arg_parser.add_argument('-o', '--save_path', required=True,
                            help='path to HDF file where model will be saved')
    args = arg_parser.parse_args()
    return args


def main():
    print(os.path.dirname(os.path.realpath(__file__)))
    logging.basicConfig(level=logging.INFO,
                       filename="%s/logs/ModelSaverLog" % os.path.dirname(os.path.realpath(__file__)))
    #log = logging.getLogger(name='GirusApp_process%d' % rank)
    log = logging.getLogger(name='ModelSaver')
    args = get_args()
    run_everything(**args.__dict__)


def create_model(num_features, layers):
    """
    Build and return a Sequential model with Dense layers given by the layers argument.
    Arguments
        num_features  (int) dimension of input
        layers     (tuple) sequence of 4-ples (units, batch_norm, activation, dropout), one per layer, such as ((64, True, 'relu', 0.5), (64, True, 'relu', 0.5), (1, False, 'sigmoid', 0.0))
    Return
        model      (Model) a compiled model
    """
    model = Sequential()
    first = True
    for (units, batch, act, dropout) in layers:
        if first is True:
            build_layer(model = model,
                        layer_type = 'Dense',
                        units = units,
                        activation = act,
                        batch_norm = batch,
                        input_dims = num_features,
                        dropout = dropout)
            first = False
        else:
            build_layer(model = model,
                        layer_type = 'Dense',
                        units = units,
                        activation = act,
                        batch_norm = batch,
                        dropout = dropout)
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    return model


def build_layer(model, layer_type, units, activation, batch_norm = True, input_dims = None, dropout = 0.0):
    if layer_type == 'Dense':
        if input_dims is None:
            model.add(Dense(units, kernel_initializer='he_uniform'))
        else:
            model.add(Dense(units, input_dim=input_dims, kernel_initializer='he_uniform'))
    if batch_norm is True:
        model.add(BatchNormalization())
    model.add(Activation(activation))
    if dropout is not 0.0:
        model.add(Dropout(dropout))


def run_everything(train_file, kmer, workers, save_path):
    log = logging.getLogger(name='ModelSaver')
    log.info(str(device_lib.list_local_devices()) + '\n')
    total_start_time = time.time()
    log.info("%s %s %d %s\n" % (train_file, kmer, workers, save_path))
    #TODO CHANGE BACK FROM RANDOM
    #train_file = "/rsgrps/bhurwitz/mattmiller899/girus_work/deep_learning/hdf_virus_3000_2_gpu.h5"
    #Step 1: walk the hdf files and get names/labels of every contig
    walk_start_time = time.time()

    #X_train, y_train, train_pos_contigs, train_neg_contigs = load_entire_hdf(train_file)
    #test_x, test_y, test_pos_contigs, test_neg_contigs = load_entire_hdf(test_file)
    train_partition, train_labels, train_pos_contigs, train_neg_contigs = walk_hdf(train_file)
    #test_partition, test_labels, test_pos_contigs, test_neg_contigs = walk_hdf(test_file)
    #log.info("X_train shape: {}\tY_train shape: {}\n".format(X_train.shape, y_train.shape))
    walk_end_time = time.time()
    log.info("Time to walk train and test: {}s\n".format(walk_end_time - walk_start_time))
    log.info("Pos neg train contigs: {} {}\n".format(train_pos_contigs, train_neg_contigs))
    #log.info("Pos neg test contigs: {} {}\n".format(test_pos_contigs, test_neg_contigs))
    #Step 2: create training generator

    training_generator = DataGenerator(train_partition, train_labels, kmer, train_file, shuffle=True)

    #Step 3: create model
    num_features = get_num_features(kmer)
    layers = ((100, True, 'relu', 0.0),
               (50, True, 'relu', 0.0),
               (25, True, 'relu', 0.0),
               (12, True, 'relu', 0.0),
               (6, True, 'relu', 0.0),
               (3, True, 'relu', 0.0),
               (1, False, 'sigmoid', 0.0))
    """

    layers = ((1024, True, 'relu', 0.0),
              (512, True, 'relu', 0.0),
              (256, True, 'relu', 0.0),
              (128, True, 'relu', 0.0),
              (64, True, 'relu', 0.0),
              (32, True, 'relu', 0.0),
              (16, True, 'relu', 0.0),
              (8, True, 'relu', 0.0),
              (4, True, 'relu', 0.0),
              (2, True, 'relu', 0.0),
              (1, False, 'sigmoid', 0.0))
    """
    #model = create_model(num_features)
    model = create_model(num_features, layers)
    #Step 4: Create callbacks
    model.summary(print_fn=lambda x: log.info(x + '\n'))
    log.info("About to start training\n")
    #, json_writer]
    #Step 5: run it
    fit_start_time = time.time()
    history = model.fit_generator(generator=training_generator, use_multiprocessing=True, workers=workers, epochs=10, verbose=2)
    #history = model.fit(X_train, y_train, batch_size=64, epochs=10, verbose=2)
    #, callbacks=callbacks)
    log.info(str(history.history.keys()))
    log.info('\n')
    log.info(str(history.history['acc']))
    log.info('\n')
    log.info(str(history.history['loss']))
    print("Saving to %s" % save_path)
    model.save(save_path)

if __name__ == "__main__":
    main()
