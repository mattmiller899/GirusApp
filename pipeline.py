import numpy as np
import argparse
from mpi4py import MPI
from Bio import SeqIO, Seq
from Bio.SeqUtils import GC
import subprocess
import logging
import pandas as pd
import itertools

def get_args():
    arg_parser = argparse.ArgumentParser()
    args = arg_parser.parse_args()
    return args


def main():
    args = get_args()
    Pipeline(**args.__dict__).run()


class Pipeline:
    def __init__(self,
                input_dir,
                k
    ):
        self.input_dir = input_dir
        self.k = k
        self.uproc_executable_fp = os.environ.get('UPROC', default='uproc')
        #TODO change to not hard-coded
        self.uproc_db_fp = '/rsgrps/bhurwitz/hurwitzlab/data/reference/uproc/pfam27ready'
        self.uproc_model_fp = '/rsgrps/bhurwitz/hurwitzlab/data/reference/uproc/model'

    def run(self):
        log = logging.getLogger(name='GirusApp')
        input_files = glob.glob(os.path.join(self.input_dir,
                                            '*.%s' % self.ext))
        out_dir = os.mkdir(os.path.join(self.input_dir,
                                            'output_dir'))
        uproc_dir = os.mkdir(os.path.join(out_dir, 'uproc_dir'))
        for input_file in input_files:
            log.info('Preprocessing file "%s"' % input_file)
            input_name = os.path.splittext(os.path.basename(input_file))[0]
            output_name = '%s.csv' % input_name
            out_fp = os.path.join(out_dir, output_name)
            count = 0
            curr_id = ''
            log.info('Running uproc')
            uproc_out = self.uproc(input_file, uproc_dir, log)
            with open(input_file, 'r') as f:
                with open(out_fp, 'w') as out:
                    for i, l in enumerate(f):
                        l = l[-1]
                        if self.ext is 'fasta':
                            if count % 2 == 1:
                                self.step_01_gc_content(l, out, curr_id)
                                self.step_02_num_orfs_and_codon_bias(l, out, i, uproc_out)
                                self.step_03_kmer_freq(l, out)
                            else:
                                curr_id = l
                        else:
                            if count % 4 == 0:
                                curr_id = l
                            elif count % 4 == 1:
                                self.step_01_gc_content(l, out, curr_id)
                                self.step_02_num_orfs_and_codon_bias(l, out, i, uproc_out)
                                self.step_03_kmer_freq(l, out)
                        count += 1
            log.info('Finished file "%s"' % input_file)
        model = self.create_model()
        feature_files = glob.glob(os.path.join(out_dir,
                                                '*.csv')
        
        for feature_file in feature_files:
            

    def step_01_gc_content(self, l, out, curr_id):
        seq = Seq(l)
        out.write('%s,%f,' % (curr_id, (GC(seq) / 100)))

 
    def step_02_num_orfs_and_codon_bias(self, l, out, i, uproc_file):
        dataset = pd.read_csv(uproc_file)
        i = i+1
        found_orfs = False
        #Get all orfs
        contig_info_dict = {}
        for _, row in dataset.iterrows():
            if row[0] == i:
                found_orfs = True
                contig_info = (int(row[3]), int(row[4]), int(row[5]))
                contig_score = float(row[7])
                if contig_info not in contig_info_dict:
                    contig_info_dict[contig_info] = contig_score
                elif contig_score > contig_info_dict[contig_info]:
                    contig_info_dict[contig_info] = contig_score
            elif row[0] > i:
                break
        if found_orfs is False:
            out.write('0,%s,' % ','.join([0.0 for x in range(64)]))
            return
        #Gt codon bias based on ORF's and sequence
        codon_counter = 0
        codon_vals = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        good_nucls = ['A', 'C', 'G', 'T']
        codon_bias_arr = [0.0 for x in range(64)]
        for key, _ in contig_info_dict.items():
            rf = key[0]
            start = key[1] - 1
            num_codons = key[2] 
            if rf < 4:
                seq = l
            else:
                seq = reverse_complement(l)
            for x in range(num_codons):
                codon = seq[start+(x*3):start+(x*3)+3]
                #If the codon has an N in it, or if the codon length is less than 3, skip this codon
                if len(codon) != 3:
                    continue
                elif codon[0] not in good_nucls or codon[1] not in good_nucls or codon[2] not in good_nucls:
                    continue
                codon_counter += 1
                codon_num = codon_vals[codon[0]] * 16 + codon_vals[codon[1]] * 4 + codon_vals[codon[2]]
                codon_bias_arr[codon_num] += 1.0
        for x in range(64):
            codon_bias_arr[x] /= codon_counter 
        out.write('%d,%s,' % (len(contig_info_dict), ','.join(codon_bias_arr)))

 
    def step_03_kmer_freq(self, l, out):
        all_kmers = gen_kmers(self.k)
        un_kmers = uniq_kmers(all_kmers)
        seq = Seq(l)
        seqUp = seq.upper()
        nkmers = len(seq) - k + 1
        kmers  = dict() 
        for i in list(range(0, nkmers - 1)):
            kmer = str(seqUp[i:i + k])

            if kmer in un_kmers:
                if kmer in kmers:
                    kmers[kmer] += 1
                else:
                    kmers[kmer] = 1
            else :
                rev = find_reverse(kmer)
                if rev in kmers:
                    kmers[rev] += 1
                else:
                    kmers[rev] = 1

        counts = [ (kmers[x]/float(nkmers)) if x in kmers else 0 for x in un_kmers ]
        out.write(",".join([str(x) for x in counts]))


    def uproc(self, input_file, uproc_dir, log):
        input_name = os.path.splittext(os.path.basename(input_file))
        output_name = '%s.csv' % input_name
        uproc_out = os.path.join(uproc_dir, output_name)
        run_cmd([
                self.uproc_executable_fp,
                '-p',
                '-o', uproc_out,
                '--pthresh', 0,
                '--othresh', 2,
                self.uproc_db_fp, 
                self.uproc_model_fp,
                input_file
                ],
                log_file=os.path.join(uproc_dir,
                                        'log')
                )
        return uproc_out


    def create_model(self):
        num_features = get_num_features(self.k)
        model = Sequential()
    
        model.add(Dense(100, input_dim=num_features, kernel_initializer='he_uniform'))
        model.add(BatchNormalization())
        model.add(Activation('relu'))
    
        model.add(Dense(50, kernel_initializer='he_uniform'))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

        model.add(Dense(25, kernel_initializer='he_uniform'))
        model.add(BatchNormalization())
        model.add(Activation('relu'))    
    
        model.add(Dense(12, kernel_initializer='he_uniform'))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

        model.add(Dense(6, kernel_initializer='he_uniform'))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

        model.add(Dense(3, kernel_initializer='he_uniform'))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

        model.add(Dense(1, kernel_initializer='he_uniform'))
        model.add(Activation('sigmoid'))
        
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        return model

def run_cmd(cmd_line_list, log_file, **kwargs):
    log = logging.getLogger(name=__name__)
    try:
        with open(log_file, 'at') as log_file:
            cmd_line_str = ' '.join((str(x) for x in cmd_line_list))
            log.info('executing "%s"', cmd_line_str)
            log_file.write('executing "{}"'.format(cmd_line_str))
            output = subprocess.run(
                cmd_line_list,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                **kwargs)
            log.info(output)
        return output
    except subprocess.CalledProcessError as c:
        logging.exception(c)
        print(c.message)
        print(c.cmd)
        print(c.output)
        raise c
    except Exception as e:
        logging.exception(e)
        print('blarg!')
        print(e)
        traceback.print_exc()
        raise e


def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement.get(base, base) for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


def find_reverse(seq):
    reverse=""
    for base in seq :
        if base == 'A':
            reverse=reverse+'T'
        elif base =='T':
            reverse=reverse+'A'
        elif base == 'C':
            reverse=reverse+'G'
        elif base =='G':
            reverse=reverse+'C'
        else :
            reverse=reverse+'X'
    return reverse


def gen_kmers(k):
    bases=['A','T','G','C']
    return [''.join(p) for p in itertools.product(bases, repeat=k)]


def uniq_kmers(kmers):
    didthat=[]
    uniq =[]
    for kmer in kmers:
        if kmer not in didthat :
            didthat.append(kmer)
            reverse=find_reverse(kmer)
            didthat.append(reverse)
            uniq.append(kmer)

    return uniq

def get_num_features(kmer):
    return int(66 + ((4 ** int(kmer)) / 2))
