import logging
import h5py
import argparse
from H5Methods import *
import random
from mpi4py import MPI
import math
import random
import time

def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-k', '--kmer', required=True,
                           help='size of kmer (can be \"no_kmer\")')
    arg_parser.add_argument('-c', '--contig', required=True, type=int,
                            help='size of contig')
    arg_parser.add_argument('-d', '--h5_path', required=True,
                            help='path to h5 output file')
    arg_parser.add_argument('-p', '--pos_dir', required=True,
                            help='path to directory containing positive (girus) feature files')
    arg_parser.add_argument('-n', '--neg_dir', required=True,
                            help='path to directory containing negative feature files')
    arg_parser.add_argument('-f', '--num_folds', required=True, type=int,
                            help='number of folds to split data into')
    args = arg_parser.parse_args()
    return args


def main():
    args = get_args()
    create_hdf(**args.__dict__)


def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def create_hdf(kmer, contig, h5_path, pos_dir, neg_dir, num_folds):
    #I do not know why, but using a flag to set the value of split_bool does not work. Need to pass a string
    create_hdf5_files_train_test(h5_path, kmer, contig, num_folds, pos_dir, neg_dir)


def create_hdf5_files_train_test(h5_path, kmer, contig, num_folds, pos_dir, neg_dir):
    start_time = time.time()
    num_processes = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    print('process {} starting up'.format(rank))
    logging.basicConfig(level=logging.INFO,
                       filename="%s/logs/HDFCreatorLog_%d" % (os.path.dirname(os.path.realpath(__file__)), rank))
    #log = logging.getLogger(name='GirusApp_process%d' % rank)
    log = logging.getLogger(name='HDF_%d' % rank)
    comm = MPI.COMM_WORLD
    #Master
    #out = open(workpath, 'w')
    #out.write('process {} starting up\n'.format(rank))
    if rank == 0:
        #out.write("Worker {} beginning master function\n".format(rank))
        #out.close()
        total_added_contigs = master(h5_path, kmer, contig, num_folds, comm, pos_dir, neg_dir)
    else:
        #out.write("Worker {} beginning worker function\n".format(rank))
        #out.close()
        worker(h5_path, num_folds, kmer, contig, comm, pos_dir, neg_dir)
    #out = open(workpath, 'a')
    #out.write("Master finished working\n")
    #out.close()
        curr_tuple = total_added_contigs[i]
        report_out.write("Iteration {}:\nPos train contigs: {}\nNeg train contigs: {}\nPos test contigs: {}\nNeg test contigs: {}\n\n\n".format(i, curr_tuple[0], curr_tuple[1], curr_tuple[2], curr_tuple[3]))
    end_time = time.time()
    log.info("Total time: {}s\n".format((end_time - start_time)))


def master(h5_path, kmer, contig, num_folds, comm, pos_dir, neg_dir):
    #out = open (outpath, 'a')
    #out.write('In master method\n')
    rank = comm.rank
    status = MPI.Status()
    tags = enum('READY', 'DONE', 'EXIT', 'START')
    num_workers = comm.size - 1
    num_pos_contigs = 0
    pos_file_tuples= []
    #Get num_pos_contigs
    #out.write("receiving pos data from workers\n")
    #out.close()
    #out = open(outpath, 'a')
    for i in range(num_workers):
        data = comm.recv(source=MPI.ANY_SOURCE, tag=tags.DONE, status=status)
        num_pos_contigs += data[0]
        pos_file_tuples.extend(data[1])
        #out.write('received\n')
    #out.close()
    #out = open(outpath, 'a')
    #Unlock workers
    #out.write('releasing workers...')
    #out.close()
    release_workers(num_workers, comm, None)
    #out = open(outpath, 'a')
    #out.write("num_pos_contigs = {}\n".format(num_pos_contigs))
    #out.write("num pos_file_names = {}\n".format(len(pos_file_tuples)))
    #out.close()
    #Get_neg_file_names
    neg_file_tuples = []
    num_neg_contigs = 0
    #out = open(utpath, 'a')
    #out.write('receiving neg data from workers...\n')
    #out.close()
    #out = open(outpath, 'a')
    for i in range(num_workers):
        data = comm.recv(source=MPI.ANY_SOURCE, tag=tags.DONE, status=status)
        #comm.Recv([data, MPI.INT], source=MPI.ANY_SOURCE, tag=tags.DONE, status=status)
        num_neg_contigs += data[0]
        neg_file_tuples.extend(data[1])
        #out.write('received\n')
    #out.write("num_neg_contigs = {}\n".format(num_neg_contigs))
    #out.write("num neg_file_tuples = {}\n".format(len(neg_file_tuples)))
    #out.close()
    #Chunk pos and neg
    pos_files_per_fold = int(len(pos_file_tuples) * (1 / num_folds))
    neg_files_per_fold = int(len(neg_file_tuples) * (1 / num_folds))
    #if not split_bool:
    random.shuffle(pos_file_tuples)
    random.shuffle(neg_file_tuples)
    print(neg_file_tuples)
    total_pos_2d_arr, pos_contig_counts_arr = chunk_total_dataset(pos_file_tuples, pos_files_per_fold, num_folds)
    total_neg_2d_arr, neg_contig_counts_arr = chunk_total_dataset(neg_file_tuples, neg_files_per_fold, num_folds)
    print("neg = " + str(total_neg_2d_arr))
    print(neg_contig_counts_arr)
    print("pos = " + str(total_pos_2d_arr))
    print(pos_contig_counts_arr)
    #MAYBE USE THIS LATER
    #write_to_split_file(split_file, total_pos_2d_arr, total_neg_2d_arr, num_folds)
    """
    else:
        total_pos_2d_arr = []
        total_neg_2d_arr = []
        pos_contig_counts_arr = []
        neg_contig_counts_arr = []
        #BIG TECHNICAL DEBT PLZ FIX FOR FUTURE MATT'S SAKE
        pos_input_dir = '/rsgrps/bhurwitz/mattmiller899/girus_work/COMPILED_FEATURES_files/girus_unsplit/{}_contigs/{}'.format(contig, kmer)
        neg_input_dir = '/rsgrps/bhurwitz/mattmiller899/girus_work/COMPILED_FEATURES_files/virus_unsplit/{}_contigs/{}'.format(contig, kmer)
        with open(split_file, 'r') as f:
            for count, l in enumerate(f):
                #Remove \n
                l = l[:-1]
                file_paths_arr = l.split('\t')
                if count % 2 == 0:
                    #curr_arr, curr_contig_count = get_contig_counts_from_tuples(file_paths_arr, pos_file_tuples, pos_input_dir)
                    curr_arr, curr_contig_count = get_contig_counts_from_tuples(file_paths_arr, pos_file_tuples, contig, kmer)
                    total_pos_2d_arr.append(curr_arr)
                    pos_contig_counts_arr.append(curr_contig_count)
                else:
                    #curr_arr, curr_contig_count = get_contig_counts_from_tuples(file_paths_arr, neg_file_tuples, neg_input_dir)
                    curr_arr, curr_contig_count = get_contig_counts_from_tuples(file_paths_arr, neg_file_tuples, contig, kmer)
                    total_neg_2d_arr.append(curr_arr)
                    neg_contig_counts_arr.append(curr_contig_count)
    """
    for i in range(num_folds):
        total_pos_str = '\t'.join(str(os.path.dirname(item)) for item in total_pos_2d_arr[i])
        total_neg_str = '\t'.join(str(os.path.dirname(item)) for item in total_neg_2d_arr[i])
        print('{}\n{}'.format(total_pos_str, total_neg_str))
    #exit()
    #out = open(outpath, 'a')
    #out.write(str(pos_contig_counts_arr))
    #out.write(str(neg_contig_counts_arr))
    #out.close()
    num_features = get_num_features(kmer)
    #Create the train/test sets, add to HDF files
    total_added_contigs = []
    # Create iterations # of HDF files, with 1 set of data in test and rest in train
    for i in range(0, num_folds):
        train_pos_arr = []
        train_neg_arr = []
        test_pos_arr = []
        test_neg_arr = []
        train_pos_contig_count = 0
        train_neg_contig_count = 0
        test_pos_contig_count = 0
        test_neg_contig_count = 0
        for j in range(0, num_folds):
            if j != i or num_folds == 1:
                train_pos_arr.extend(total_pos_2d_arr[j])
                train_neg_arr.extend(total_neg_2d_arr[j])
                train_pos_contig_count += pos_contig_counts_arr[j]
                train_neg_contig_count += neg_contig_counts_arr[j]
            else:
                test_pos_arr.extend(total_pos_2d_arr[j])
                test_neg_arr.extend(total_neg_2d_arr[j])
                test_pos_contig_count += pos_contig_counts_arr[j]
                test_neg_contig_count += neg_contig_counts_arr[j]
        #out = open(outpath, 'a')
        #out.write('\ntrain_pos_contig_count = {}\ntrain_neg_contig_count = {}\ntest_pos_contig_count = {}\ntest_neg_contig_count = {}\n'.format(train_pos_contig_count, train_neg_contig_count, test_pos_contig_count, test_neg_contig_count))
        #out.close()
        send_tuple = (train_pos_arr, train_neg_arr, test_pos_arr, test_neg_arr, i, h5_path, num_features, train_pos_contig_count,
                        train_neg_contig_count, test_pos_contig_count, test_neg_contig_count)
        release_workers(num_workers, comm, send_tuple)
        added_contigs = make_hdf_files(train_pos_arr, train_neg_arr, test_pos_arr, test_neg_arr, train_pos_contig_count,
                                       train_neg_contig_count, test_pos_contig_count, test_neg_contig_count, i, h5_path,
                                       num_features, rank, comm)
        added_train_pos = 0
        added_train_neg = 0
        added_test_pos = 0
        added_test_neg = 0
        data = None
        #out = open(outpath, 'a')
        #out.write('receiving hdf stats from workers...\n')
        #out.close() 
        for j in range(num_workers):
            data = comm.recv(source=MPI.ANY_SOURCE, tag=tags.DONE, status=status)
            added_train_pos += data[0]
            added_train_neg += data[1]
            added_test_pos += data[2]
            added_test_neg += data[3]
        #out = open(outpath, 'a')
        #out.write('checking if data from workers is good...\n')
        if added_train_pos != train_pos_contig_count:
            print("ERROR WITH TRAIN_POS!!!")
            #out.write('train_pos not correct!!!\n')
        if added_train_neg != train_neg_contig_count:
            #out.write('train_neg not correct!!!\n')
            print("ERROR WITH TRAIN_NEG!!!")
        if added_test_pos != test_pos_contig_count:
            #out.write('test_pos not correct!!!\n')
            print("ERROR WITH TEST_POS!!!")
        if added_test_neg != test_neg_contig_count:
            #out.write('test_neg not correct!!!\n')
            print("ERROR WITH TEST_NEG!!!")
        #out.write('File {} finished\n'.format(i))
        #out.close()
        added_tuple = (added_train_pos, added_train_neg, added_test_pos, added_test_neg)
        total_added_contigs.append(added_tuple)
    return total_added_contigs
 

def release_workers(num_workers, comm, send_tup):
    tags = enum('READY', 'DONE', 'EXIT', 'START')
    #out = open(outpath, 'a')
    #out.write('in release workers\n')
    #out.close()
    #out = open(outpath, 'a')    
    for i in range(num_workers):
        #out.write("releasing worker {}\n".format(i + 1))
        comm.isend(send_tup, dest=i + 1, tag=tags.READY)
    #out.close()

"""
#DEPRECIATED: MAY UYSE LATER
def write_to_split_file(split_file, total_pos_2d_arr, total_neg_2d_arr, num_folds):
    with open(split_file, 'w') as f:
        for i in range(num_folds):
            #total_pos_str = '\t'.join(str(os.path.basename(os.path.dirname(item))) for item in total_pos_2d_arr[i])
            #total_neg_str = '\t'.join(str(os.path.basename(os.path.dirname(item))) for item in total_neg_2d_arr[i])
            total_pos_str = '\t'.join(str(os.path.dirname(item)) for item in total_pos_2d_arr[i])
            total_neg_str = '\t'.join(str(os.path.dirname(item)) for item in total_neg_2d_arr[i])
            f.write('{}\n{}\n'.format(total_pos_str, total_neg_str))


def get_contig_counts_from_tuples(file_paths_arr, file_tuples, input_dir):
    curr_arr = []
    curr_contig_count = 0
    for tmp_file in file_paths_arr:
        tmp_file = '{}/{}'.format(input_dir, tmp_file)
        for tmp_tuple in file_tuples:
            if tmp_file in tmp_tuple:
                curr_arr.append('{}/c0-{}'.format(tmp_file, tmp_tuple[1]))
                curr_contig_count += tmp_tuple[1]
                break
    return curr_arr, curr_contig_count


def get_contig_counts_from_tuples(file_paths_arr, file_tuples, contig, kmer):
    curr_arr = []
    curr_contig_count = 0
    old_contig = os.path.basename(file_paths_arr[0].split("_contigs")[0])
    old_kmer = os.path.dirname(file_paths_arr[0].split("_contigs/")[1])
    for tmp_file in file_paths_arr:
        tmp_file = tmp_file.replace("{}_contigs/{}".format(old_contig, old_kmer), "{}_contigs/{}".format(contig, kmer))
        for tmp_tuple in file_tuples:
            if tmp_file in tmp_tuple:
                curr_arr.append('{}/c0-{}'.format(tmp_file, tmp_tuple[1]))
                curr_contig_count += tmp_tuple[1]
                break
    return curr_arr, curr_contig_count
"""

def worker(h5_path, num_folds, kmer, contig, comm, pos_dir, neg_dir):
    #out = open(outpath, 'a')
    status = MPI.Status()
    tags = enum('READY', 'DONE', 'EXIT', 'START')
    #Subtract 1 for the master
    rank = comm.rank - 1
    num_processes = comm.size - 1
    #out.write("in worker method, about to get pos data\n")
    #out.close()
    pos_glob = "%s/*" % (pos_dir)
    num_pos_contigs, pos_file_tuples = get_pos_contig_counts(pos_glob, num_processes, rank)
    send_tup = (num_pos_contigs, pos_file_tuples)
    comm.isend(send_tup, dest=0, tag=tags.DONE)
    #out = open(outpath, 'a')
    #Wait for master/other workers to move forward
    #out.write("Worker {} waiting after sending pos data...\n".format(rank))
    #out.close()
    data = comm.recv(source=0, tag=tags.READY)
    #comm.Recv([data, MPI.INT], source=0, tag=tags.READY)

    #out = open(outpath, 'a')
    #out.write("Worker {} released\n".format(rank))
    #out.close()
    #Get_neg_file_names
    num_neg_contigs, neg_file_tuples = get_neg_file_names(neg_dir, contig, kmer, num_processes, rank)
    send_tup = (num_neg_contigs, neg_file_tuples)
    comm.isend(send_tup, dest=0, tag=tags.DONE)
    #exit()
    #out = open(outpath, 'a')
    #out.write("Sent the neg data\n".format(rank))
    #out.close()
    for j in range(num_folds):
        #out = open(outpath, 'a')
        #Waiting to start adding shit to the HDF files
        #out.write("Worker {} waiting...\n".format(rank))
        #out.close()
        (train_pos_arr, train_neg_arr, test_pos_arr, test_neg_arr, i, h5_path, num_features, train_pos_contig_count,
         train_neg_contig_count, test_pos_contig_count, test_neg_contig_count) = comm.recv(source=0, tag=tags.READY)
        #out = open(outpath, 'a')
        #out.write('Worker unlocked with hdf data\n')
        #out.close()
        make_hdf_files(train_pos_arr, train_neg_arr, test_pos_arr, test_neg_arr, train_pos_contig_count,
                        train_neg_contig_count, test_pos_contig_count, test_neg_contig_count, i, h5_path,
                        num_features, rank, comm)
        #out = open(outpath, 'a')
        #out.write('Worker finished making {} hdf file\n'.format(j))
        #out.close()
    #out = open(outpath, 'a')
    #out.write("Worker totally finally finished")
    #out.close()
    exit()


#Opens HDF files, passes to add_pos_neg_datasets_to_file
def make_hdf_files(train_pos_arr, train_neg_arr, test_pos_arr, test_neg_arr, train_pos_contig_count,
                   train_neg_contig_count, test_pos_contig_count, test_neg_contig_count, iteration, h5_path,
                   num_features, worker_rank, comm):
    tags = enum('READY', 'DONE', 'EXIT', 'START')
    #out = open(outpath, 'a')
    #out.write('starting make_hdf_files\n')
    h5_name = os.path.splitext(h5_path)[0]
    train_hdf = "{}_{}_train.h5".format(h5_name, iteration)
    test_hdf = "{}_{}_test.h5".format(h5_name, iteration)
    train_out = h5py.File(train_hdf, "w", driver='mpio', comm=MPI.COMM_WORLD)
    test_out = h5py.File(test_hdf, "w", driver='mpio', comm=MPI.COMM_WORLD)
    pos_train_dset = train_out.create_dataset(name = 'Girus',
                                             shape = (train_pos_contig_count, num_features),
                                             dtype=np.float64)
    neg_train_dset = train_out.create_dataset(name = 'Not Girus',
                                             shape = (train_neg_contig_count, num_features),
                                             dtype=np.float64)
    pos_test_dset = test_out.create_dataset(name = 'Girus',
                                             shape = (test_pos_contig_count, num_features),
                                             dtype=np.float64)
    neg_test_dset = test_out.create_dataset(name = 'Not Girus',
                                             shape = (test_neg_contig_count, num_features),
                                             dtype=np.float64)

    #Note: Worker rank is for workers, (0 through num_processes) - 1. Master has rank 0, but worker 0 has worker rank 0, real rank 1
    if MPI.COMM_WORLD.rank == 0:
        #out.write('master returning from make_hdf_files\n')
        #out.close()
        train_out.close()
        test_out.close()
        return
    #out.write('datasets created, adding to datasets\n')
    #out.write('starting train_pos\n')
    #out.close()
    train_pos_contigs = add_to_datasets(train_pos_arr, pos_train_dset, worker_rank)
    #out = open(outpath, 'a')
    #out.write('starting train_neg\n')
    #out.close()
    train_neg_contigs = add_to_datasets(train_neg_arr, neg_train_dset, worker_rank)
    #out = open(outpath, 'a')
    #out.write('starting test_pos\n')
    #out.close()
    test_pos_contigs = add_to_datasets(test_pos_arr, pos_test_dset, worker_rank)
    #out = open(outpath, 'a')
    #out.write('starting test_neg\n')
    #out.close()
    test_neg_contigs = add_to_datasets(test_neg_arr, neg_test_dset, worker_rank)
    train_out.close()
    test_out.close()
    send_tup = (train_pos_contigs, train_neg_contigs, test_pos_contigs, test_neg_contigs)
    #out = open(outpath, 'a')
    #out.write('Finished, sending data to master\n')
    #out.close()
    comm.isend(send_tup, dest=0, tag=tags.DONE)
    return


def add_to_datasets(data, dset, rank):
    #out = open(outpath, 'a')
    #out.write('adding samples to datasets\n')
    #out.close()
    num_contigs_added = 0
    num_processes = MPI.COMM_WORLD.size - 1
    contig_count = 0
    #For each path/c1-100, get path and end_contig and open a file
    for i, entry in enumerate(data):
        file_path = os.path.dirname(entry)
        contigs = os.path.basename(entry)
        end_contig = int(contigs.split('c')[-1].split('-')[1])
        if i % num_processes != rank:
            contig_count += end_contig
            continue
        with open(file_path, "r") as infile:
            for line_count, l in enumerate(infile):
                if line_count == end_contig:
                    break
                larr = np.fromstring(l, sep=",")
                dset[contig_count + line_count] = larr
                num_contigs_added += 1
        contig_count += end_contig
    return num_contigs_added


#DEPRECIATED: Creates pos neg groups, passes to add_datasets_to_group
def add_pos_neg_datasets_to_file(pos, neg, out, num_features):
    pos_group = out.create_group("Girus")
    neg_group = out.create_group("Not Girus")
    pos_contigs = add_datasets_to_pos_group(pos_group, pos, num_features)
    neg_contigs = add_datasets_to_neg_group(neg_group, neg, num_features)
    return pos_contigs, neg_contigs

        
#DEPRECIATEDCreates org groups, puts datasets in them, returns # of contigs added
def add_datasets_to_pos_group(group, data, num_features):
    num_contigs = 0
    for entry in data:
        file_path = os.path.dirname(entry)
        contigs = os.path.basename(entry)
        end_contig = contigs.split('c')[-1].split('-')[1]
        org_name = os.path.basename(file_path).split(".")[0]
        org_group = group.create_group(org_name)
        with open(file_path, "r") as infile:
            for line_count, l in enumerate(infile):
                if line_count == int(end_contig):
                    break
                larr = np.fromstring(l, sep = ",")
                contig_dset = org_group.create_dataset(name = "c%d" % line_count,
                                                       data = larr,
                                                       shape=(1, num_features),
                                                       dtype=np.float64,
                                                       compression='gzip')
                num_contigs += 1
    return num_contigs


#DEPRECIATED: Creates org groups, puts datasets in them, returns # of contigs added
def add_datasets_to_neg_group(group, data, num_features):
    num_contigs = 0
    for file_path in data:
        contigs = reservoir_sampling(file_path, 100)
        org_name = os.path.basename(file_path).split(".")[0]
        org_group = group.create_group(org_name)
        for line_count, l in enumerate(contigs):
            larr = np.fromstring(l, sep = ",")
            contig_dset = org_group.create_dataset(name = "c%d" % line_count,
                                                   data = larr,
                                                   shape=(1, num_features),
                                                   dtype=np.float64,
                                                   compression='gzip')
            num_contigs += 1
    return num_contigs

#DEPRECIATED: Randomly selects num_contigs from a file_path of unknown length
def reservoir_sampling(file_path, num_contigs):
    sample = []
    with open(file_path, 'r') as fp:
        for i,line in enumerate(fp):
            if i < num_contigs:
                sample.append(line)
            elif i >= num_contigs and random.random() < num_contigs/float(i+1):
                replace = random.randint(0,len(sample)-1)
                sample[replace] = line
    return sample


#Takes list of file names and # of contigs needed per fold, return 2d array contain files in each fold + contig# at the end and array containing contig counts for each fold
def chunk_total_dataset(file_tuples, num_files_per_fold, num_folds):
    ret_arr = []
    count_arr = []
    curr_fold = []
    #TODO CHANGE BACK
    fold_count = 0
    contig_count = 0
    for i, (input_file, curr_num_contigs) in enumerate(file_tuples):
        if i % num_files_per_fold == 0 and i != 0 and fold_count != num_folds-1:
            ret_arr.append(curr_fold)
            count_arr.append(contig_count)
            curr_fold = []
            fold_count += 1
            contig_count = 0
        curr_fold.append("{}/c{}-{}".format(input_file, 0, curr_num_contigs))
        contig_count += curr_num_contigs
    #Handle last set
    ret_arr.append(curr_fold)
    count_arr.append(contig_count)
    return ret_arr, count_arr


#DEPRECIATED: Takes list of file names, return 2d array contain files in each fold
def chunk_total_dataset_neg(file_names, iterations):
    ret_arr = []
    curr_fold = []
    #TODO CHANGE BACK
    num_folds = 0
    num_files = len(file_names)
    files_per_fold = int(num_files / iterations)
    for input_count, input_file in enumerate(file_names):
        if num_folds == 5:
            return ret_arr
        curr_fold.append("{}".format(input_file))
        if input_count % files_per_fold == 0 and input_count != 0:
            ret_arr.append(curr_fold)
            curr_fold = []
            num_folds += 1
    #Handle last set
    ret_arr.append(curr_fold)
    return ret_arr


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def get_num_features(kmer):
    if kmer == "no_kmer":
        return int(66)
    else:
        return int(66 + ((4 ** int(kmer)) / 2))


def get_pos_contig_counts(input_glob, num_processes, rank):
    file_tuples = []
    num_contigs = 0
    for i, input_file in enumerate(sorted(glob.glob(input_glob))):
        if os.path.isdir(input_file) or i % num_processes != rank:
            continue
        tmp_file_len = file_len(input_file)
        num_contigs += tmp_file_len
        tmp_tuple = (input_file, tmp_file_len)
        file_tuples.append(tmp_tuple)
    return num_contigs, file_tuples


def get_neg_file_names(neg_dir, contig, kmer, num_processes, rank):
    file_tuples = []
    num_contigs = 0
    neg_glob = "%s/*" % (neg_dir)
    for i, input_file in enumerate(glob.glob(neg_glob)):
        if os.path.isdir(input_file) or i % num_processes != rank:
            continue
        tmp_file_len = file_len(input_file)
        num_contigs += tmp_file_len
        tmp_tuple = (input_file, tmp_file_len)
        file_tuples.append(tmp_tuple)
    return num_contigs, file_tuples

if __name__ == "__main__":
    main()

