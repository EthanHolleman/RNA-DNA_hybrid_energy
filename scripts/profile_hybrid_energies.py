import numpy as np
from Bio import SeqIO
from statistics import mean
import csv


ENERGY_DICT = {
    'CC': -0.36,
    'CG': -0.16,
    'CT': -0.1,
    'CA': -0.06,
    'GC': 0.97,
    'GG': 0.34,
    'GT': 0.45,
    'GA': 0.38,
    'TC': -0.12,
    'TG': -0.16,
    'TT': 0.6,
    'TA': -0.12,
    'AC': .45,
    'AG': .5,
    'AT': .28,
    'AA': .8
}

AVERAGE_ENERGY = mean(ENERGY_DICT.values())


def read_fasta_entries(fasta_path):
    return SeqIO.parse(fasta_path, 'fasta')


def hybrid_energy(seq):
    di_nucs = zip(seq[:-1], seq[1:])
    vals = np.full(len(seq) - 1, 0.0)
    for i, di_nuc in enumerate(di_nucs):
        try:
            vals[i] = ENERGY_DICT[''.join(di_nuc).upper()]
        except KeyError as e:
            vals[i] = AVERAGE_ENERGY
    return vals


def windowed_hybrid_energy(record, window_size, dist_between_windows, chr_name, output_file):
    window_start_position = 0
    window_end_position = window_size
    seq_len = len(record)
    previous_window_vals = None
    final_window = False  # detect when need to break loop
    with open(output_file, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        while True:
            current_window_seq = str(record.seq[window_start_position:
                                                window_end_position]
                                    )
            
            current_window_vals = seq_hybrid_energy(
                current_window_seq, previous_window_vals,
                dist_between_windows
            )
            print(len(current_window_vals))
        
            previous_window_vals = current_window_vals

            writer.writerow(
                (chr_name, window_start_position,
                window_start_position+dist_between_windows,
                round(current_window_vals.sum(), 3)
                )
            )
            # set up for the next window
            window_start_position += dist_between_windows
            window_end_position += dist_between_windows
            if window_end_position >= seq_len:
                if final_window:
                    break
                else:
                    window_end_position = seq_len-1
                    final_window = True


def seq_hybrid_energy(current_seq, prev_seq_vals, dist_between_windows):
    if type(prev_seq_vals) != np.ndarray:
        return hybrid_energy(current_seq)
    else:
        overlap = prev_seq_vals[dist_between_windows:]
        # get current_seq that is not overlapping with prev_seq_vals
        non_overlap_seq = current_seq[len(overlap):]
        non_overlap_vals = hybrid_energy(non_overlap_seq)
        return np.concatenate((overlap, non_overlap_vals))


def hybrid_energy_dict_to_bed_graph(hed):
    with open('test.bed', 'w') as handle:
        import csv
        writer = csv.writer(handle, delimiter='\t')
        for key, value in hed.items():
            writer.writerow((key[0], key[1], value))


def main():
    # Should be run as a snakemake script
    # chr_name = str(snakemake.params.chr_name)
    # input_fasta = str(snakemake.input)
    # output_bed = str(snakemake.output)
    chr_name = 'chr1'
    input_fasta='/home/ethan/Documents/github/hybrid_energy/rawdata/hg19/chr2.fa'
    output_bed = 'chr10.hybrid_init.bedgraph'

    # really should be only one entry
    records = read_fasta_entries(input_fasta)
    for record in records:
        windowed_hybrid_energy(
            record, 50,
            8, chr_name,
            output_bed)


if __name__ == '__main__':
    #main()
    print(sum(hybrid_energy('N'*50)))
