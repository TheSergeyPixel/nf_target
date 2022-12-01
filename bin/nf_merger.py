import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='nf metrics merger')

parser.add_argument('-i', '--input', required=True, help='Input as path to directory with metrics files', type=str)


args = parser.parse_args()

target_dir = os.path.abspath(args.input)


def get_alignment_metrics(sample_name):
    reads_value = int(pd.read_csv(f'{target_dir}/{sample_name}_alignment_summary_metrics.txt',
                                  sep='\t', comment='#', nrows=3)['TOTAL_READS'][2])

    nucleotides_value = int(pd.read_csv(f'{target_dir}/{sample_name}_alignment_summary_metrics.txt',
                                        sep='\t', comment='#', nrows=3)['PF_ALIGNED_BASES'][2])

    read_size_value = pd.read_csv(f'{target_dir}/{sample_name}_alignment_summary_metrics.txt',
                                  sep='\t', comment='#', nrows=3)['MEAN_READ_LENGTH'][2]

    return reads_value, nucleotides_value, read_size_value


def get_hs_metrics(sample_name):
    depth_value = float(pd.read_csv(f'{target_dir}/{sample_name}_HS.txt',
                                    sep='\t', comment='#', nrows=1)['MEAN_TARGET_COVERAGE'][0])

    target_enrichment_value = float(pd.read_csv(f'{target_dir}/{sample_name}_HS.txt',
                                                sep='\t', comment='#', nrows=1)['PCT_SELECTED_BASES'][0] * 100)

    return depth_value, target_enrichment_value


def get_duplicates_metrics(sample_name):
    duplicates_value = float(pd.read_csv(f'{target_dir}/{sample_name}_dup_metrics.txt', nrows=2,
                                         sep='\t', comment='#')['PERCENT_DUPLICATION'][0] * 100)

    return duplicates_value


def get_vcf_stats(sample_name):
    variants_value = pd.read_csv(f'{target_dir}/{sample_name}_vcf_stats.txt', sep=':', comment='#', header=None).T
    variants_value.columns = variants_value.iloc[0].apply(lambda x: x.strip('               '))
    variants_value = variants_value[1:]
    variants_value = int(variants_value['Passed Filters'])

    return variants_value


def get_evenness(sample_name, depth):
    evenness_df = pd.read_csv(f'{target_dir}/{sample_name}_per_target.txt', sep='\t')
    evenness_value = round(len(evenness_df[evenness_df['mean_coverage'] > 0.2 * depth]) / len(evenness_df) * 100, 1)

    return evenness_value


def main():
    df = pd.DataFrame()

    sample_names = [a.strip('.pdf') for a in list(os.listdir(target_dir)) if a.endswith('.pdf')]

    sample_list = []
    reads_list = []
    nucleotides_list = []
    depth_list = []
    variants_list = []
    duplicates_list = []
    mean_size_list = []
    target_enrichment_list = []
    evenness_list = []

    for sample in sample_names:
        reads, nucleotides, read_size = get_alignment_metrics(sample)

        depth, target_enrichment = get_hs_metrics(sample)

        duplicates = get_duplicates_metrics(sample)

        variants = get_vcf_stats(sample)

        evenness = get_evenness(sample, depth)

        sample_list.append(sample)
        reads_list.append(reads)
        nucleotides_list.append(nucleotides)
        depth_list.append(depth)
        variants_list.append(variants)
        duplicates_list.append(duplicates)
        mean_size_list.append(read_size)
        target_enrichment_list.append(target_enrichment)
        evenness_list.append(evenness)

    df['Sample'] = sample_list
    df['Total_Reads'] = reads_list
    df['Total_Nucleotides'] = nucleotides_list
    df['Depth'] = depth_list
    df['Total_Variants'] = variants_list
    df['Duplicates_%'] = duplicates_list
    df['Mean_Read_Size'] = mean_size_list
    df['Target_Read_Enrichment_%'] = target_enrichment_list
    df['Evenness_Of_Coverage'] = evenness_list

    return df.to_csv('Metrics_results.tsv', sep='\t', index=None)


if __name__ == "__main__":
    main()
