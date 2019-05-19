from matplotlib.patches import Rectangle  # Rectangle is used despite it being greyed out in pycharm
from matplotlib.collections import PatchCollection
# https://stackoverflow.com/questions/21784641/installation-issue-with-matplotlib-python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
from collections import defaultdict
import sys
class InvestigateBIDSeqs():
    def __init__(self):
        self.base_data_dir = '/Users/humebc/Google_Drive/projects/tara/testing_BID_adapters/data_testing'
        self.seq_absolute_output_path = '/Users/humebc/Google_Drive/projects/tara/testing_BID_adapters/data_loading_outputs'
        self.seqs = ['C1', 'C1d', 'C1b', 'C42a', 'C42g', 'C42.2']
        self.sample_names = ['BUR_AAYN', 'BUR_AEZY', 'BUR_AFAE', 'BUM_ABBZ', 'BUM_AAFA', 'BUM_ABBW']
        self.df = pd.DataFrame(index=self.sample_names, columns=self.seqs)
        self.seqs_dict = {
            'C1'   :'AATGGCCTCCTGAACGTGCGTTGCACTCTTGGGATTTCCTGAGAGTATGTCTGCTTCAGTGCTTAACTTGCCCCAACTTTGCAAGCAGGATGTGTTTCTGCCTTGCGTTCTTATGAGCTATTGCCCTCTGAGCCAATGGCTTGTTAATTGCTTGGTTCTTGCAAAATGCTTTGCGCGCTGTTATTCAGGTTTCTACCTTCGTGGTTTTACTTGAGTGACGCTGCTCATGCTTGCAACCGCTGGGATGCAGGTGCATGCCTCTA',
            'C1d'  :'AATGGCCTCCTGAACGTGCGTTGCACTCTTGGGATTTCCTGAGAGTATGTCTGCTTCAGTGCTTAACTTGCCCCAACTTTGCAAGCAGGATGTGTTTCTGCCTTGCGTTCTTATGAGCTATTGCCCTCTGAGCCAATGGCTTGTGAATTGCTTGGTTCTTGCAAAATGCTTTGCGCGCGGTTATTCAGGTTTCTACCTTCGTGGTTTTACTTGAGTGACGCTGCTCATGCTTGCAACCGCTGGGATGCAGGTGCATGCCTCTA',
            'C1b'  :'AATGGCCTCCTGAACGTGCGTTGCACTCTTGGGATTTCCTGAGAGTATGTCTGCTTCAGTGCTTAACTTGCCCCAACTTTGCAAGCAGATGTGTTTCTGCCTTGCGTTCTTATGAGCTATTGCCCTCTGAGCCAATGGCTTGTTAATTGCTTGGTTCTTGCAAAATGCTTTGCGCGCTGTTATTCAGGTTTCTACCTTCGTGGTTTTACTTGAGTGACGCTGCTCATGCTTGCAACCGCTGGGATGCAGGTGCATGCCTCTA',
            'C42a' :'AATGGCCTCCTGAACGTGCGTTGCACTCTTGGGATTTCCTGAGAGTATGTCTGCTTCAGTGCTTAACTTGCCCCAACTTTGCAAGCAGGATGTGTTTCTGCCTTGCGTTCTTATGCGCTATTGCCCTCTGAGCCAATGGCTTGTGAATTGCTTGGTTCTTGCAAAATGCTTTGCGCGCTGTTATTCAGGTTTCTACCTTCGTGGTTTTACTTGAGTGACGCTGCTCATGCTTGCAACCGCTGGGATGCAGGTGCATGCCTCTA',
            'C42g' :'AATGGCCTCCTGAACGTGCGTTGCACTCTTGGGATTTCCTGAGAGTATGTCTGCTTCAGTGCTTAACTTGCCCCAACTTTGCAAGCAGGACGTGTTTCTGCCTTGCGTTCTTATGCGCTATTGCCCTCTGAGCCAATGGCTTGTGAATTGCTTGGTTCTTGCAAAATGCTTTGCGCGCTGTTATTCAGGTTTCTACCTTCGTGGTTTTACTTGAGTGACGCTGCTCATGCTTGCAACCGCTGGGCTGCAGGTGCATGCCTCTA',
            'C42.2':'AATGGCCTCCTGAACGTGCGTTGCACTCTTGGGATTTCCTGAGAGTATGTCTGCTTCAGTGCTTAACTTGCCCCAACTTTGCAAGCAGGATGTGTTTCTGCCTTGCGTTCTTATGAGCTATTGCCCTCTGAGCCAATGGCTTGTGAATTGCTTGGTTCTTGCAAAATGCTTTGCGCGCTGTTATTCAGGTTTCTACCTTCGTGGTTTTACTTGAGTGACGCTGCTCATGCTTGCAACCGCTGGGATGCAGGTGCATGCCTCTA'
        }
        self.meta_info_path = '/Users/humebc/Google_Drive/projects/tara/testing_BID_adapters/bid_meta_info.csv'

    def investigate(self):
        for sample in self.sample_names:
            print(f'Proceessing {sample}')
            sample_dir = os.path.join(self.base_data_dir, sample)
            fasta_path = os.path.join(sample_dir, 'filetrim.contigs.fasta')
            fasta_as_list = None
            with open(fasta_path, 'r') as f:
                fasta_as_list = [line.rstrip().lstrip() for line in f]
            num_seqs = len(fasta_as_list)/2
            dd = defaultdict(int)
            for i in range(1, len(fasta_as_list), 2):
                count = 0
                for k, v in self.seqs_dict.items():
                    if v in fasta_as_list[i]:
                        dd[k] += 1
                        count +=1
                if count > 1:
                    sys.exit('seq contined more than one of the seqs')
            for seq in self.seqs:
                if seq in dd:
                    self.df.at[sample, seq] = dd[seq]/num_seqs
                else:
                    self.df.at[sample, seq] = 0


ibids = InvestigateBIDSeqs()
