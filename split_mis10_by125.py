import os
import sys
import pysam


def split_sam_files_by_insert_size(sample_path, sample_name_list_path,
                                   output_path):
    with open(sample_name_list_path, 'r') as f:
        sample_names = [line.strip() for line in f if line.strip()]

    total_samples = len(sample_names)
    processed_samples = 0

    for sample_name in sample_names:
        bam_file_path = os.path.join(sample_path, sample_name,
                                     f"{sample_name}.mis.10.bam")
        if os.path.getsize(bam_file_path) > 1:
            print(f'Processing file {bam_file_path}')

            sam_file_path_4 = os.path.join(output_path,
                                           f"{sample_name}.mis.10.S125.sam")
            sam_file_path_5 = os.path.join(output_path,
                                           f"{sample_name}.mis.10.L125.sam")

            with pysam.AlignmentFile(bam_file_path, 'rb') as sam_ali, \
                 pysam.AlignmentFile(sam_file_path_4, 'w', template=sam_ali) as sam_4, \
                 pysam.AlignmentFile(sam_file_path_5, 'w', template=sam_ali) as sam_5:

                for read in sam_ali.fetch():
                    insert_size = abs(int(read.template_length))
                    if insert_size < 125:
                        sam_4.write(read)
                    elif insert_size >= 125:
                        sam_5.write(read)

            processed_samples += 1
            print(
                f'{processed_samples}/{total_samples} samples have been processed'
            )
            if processed_samples == total_samples:
                print('All samples have been processed')
        else:
            print(f'File {bam_file_path} is too small or does not exist')


def main():
    sample_path = sys.argv[1]
    sample_name_list_path = os.path.join(sample_path,'sample_name.txt')
    output_path = os.path.join(sample_path,
                               'fragment_study/mis_10_splitedBy125',
                               '125_splited_sam_file')
    os.makedirs(output_path, exist_ok=True)

    split_sam_files_by_insert_size(sample_path, sample_name_list_path,
                                   output_path)


if __name__ == '__main__':
    main()
