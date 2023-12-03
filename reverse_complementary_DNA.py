import sys

def complementary(seq):
    base_dic = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N','D':'D','I':'I'}
    seq = seq[::-1]
    comp_seq = ''
    for base in seq:
        comp_seq += base_dic[base]
    return comp_seq

if __name__ == '__main__':
    seq = sys.argv[1]
    print(complementary(seq))