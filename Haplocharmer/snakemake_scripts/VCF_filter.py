# VCF_filter

# Import modules
import sys
import os
import gzip
import datetime
import argparse
import pandas as pd
from collections import defaultdict


# Function main
def main():
    # Create the parser object
    parser = argparse.ArgumentParser(description='VCF_filter.py is a python script to apply filters on a VCF '
                                                 'resulting from the merging of gVCF files')
    # Add arguments to the parser
    parser.add_argument('-i', '--input', type=str, help='Input VCF file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output VCF file prefix', required=True)
    parser.add_argument('--compress', action='store_true', help='Compress output using gzip')
    parser.add_argument('--min_gt_depth', type=int, default=10,
                        help='Minimum read depth (DP) of the genotype datapoint')
    parser.add_argument('--max_gt_depth', type=int, default=1000,
                        help='Maximum read depth (DP) of the genotype datapoint')
    parser.add_argument('--min_hap_depth', type=int, default=3,
                        help='Minimum haplotype depth (AD) of the genotype datapoint over all individuals')
    parser.add_argument('--min_hap_freq', type=float, default=0.04,
                        help='Minimum haplotype frequency  (AD/DP) of the genotype datapoint over all individuals')
    # Parse the arguments
    args = parser.parse_args()
    # Apply checks
    output_name, input_read_mode, input_read_pd_mode, output_suffix, output_write_mode, \
        output_write_pd_mode = process_checks(args)
    # Open VCF file to write the header, get column names and phasing set sizes
    with (gzip.open if args.input.endswith(".vcf.gz") else open)(args.input, input_read_mode) as input_file, (
            gzip.open if args.compress else open)(output_name + output_suffix, output_write_mode) as output_file:
        # Generate header and column names
        file_colnames = generate_header(input_file, output_file)
        # Create a dictionary with the dimension of each phasing set
        ps_sizes = get_ps_sizes(input_file)
    # Filter the VCF phasing set by phasing set
    process_filters_by_ps(args, ps_sizes, file_colnames, output_name, input_read_pd_mode, output_suffix,
                          output_write_pd_mode)


# Function to filter the VCF phasing set by phasing set
def process_filters_by_ps(args, ps_sizes, vcf_colnames, output_name, vcf_read_pd_mode, vcf_suffix, vcf_write_pd_mode):
    # Open VCF file
    vcf_reader = pd.read_table(args.input, sep="\t", iterator=True, header=None, comment='#',
                               compression=vcf_read_pd_mode)
    # Iterate through phasing sets
    ps_incr = 1
    for ps_key, ps_value in ps_sizes.items():
        # Get chunk of VCF corresponding to the phasing set
        ps_vcf = (vcf_reader.get_chunk(ps_value + ps_incr))
        ps_incr = 0
        # Set column names
        ps_vcf.columns = vcf_colnames
        # Apply filtering of haplotypes
        ps_vcf_filt = filter_hap(ps_vcf, args.min_gt_depth, args.max_gt_depth, args.min_hap_depth,
                                 args.min_hap_freq)
        # Apply filtering of non-informative variants and reset ALT alleles
        ps_vcf_filt_curated = ps_vcf_filt.apply(filter_reset_variant, axis=1, result_type='broadcast').dropna()
        # Write filtered vcf block
        ps_vcf_filt_curated.to_csv(output_name + vcf_suffix, sep="\t", index=False, mode="a", header=False,
                                   compression=vcf_write_pd_mode)


# Function to apply checks
def process_checks(args):
    # Test for the presence of the VCF file
    if not os.path.isfile(args.input):
        sys.exit("There is no VCF file at this path")
    # Test for VCF file type
    if args.input.endswith(".vcf.gz"):
        input_read_mode = "rt"
        input_read_pd_mode = "gzip"
    elif args.input.endswith(".vcf"):
        input_read_mode = "r"
        input_read_pd_mode = None
    else:
        sys.exit("Invalid input type, please use .vcf or .vcf.gz files")
    # Create output directory if it does not exist already
    output_dir = os.path.dirname(args.output)
    if len(output_dir) > 0:
        os.makedirs(output_dir, exist_ok=True)
    # Remove .vcf, .hpa ang .gz extensions from output if any
    output_name, output_extension = os.path.splitext(args.output)
    while output_extension == ".vcf" or output_extension == ".hpa" or output_extension == ".txt" or \
            output_extension == ".gz":
        output_name, output_extension = os.path.splitext(output_name)
    # Build VCF suffix
    output_suffix = ".vcf"
    # Test if gVCF compression is needed
    if args.compress:
        output_suffix += ".gz"
        output_write_mode = "wt"
        output_write_pd_mode = "gzip"
    else:
        output_write_mode = "w"
        output_write_pd_mode = None
    # Return
    return output_name, input_read_mode, input_read_pd_mode, output_suffix, output_write_mode, output_write_pd_mode


# Function to generate the header
def generate_header(input_file, output_file):
    # Create header list
    file_header = ['##fileformat=VCFv4.2']
    # Get the current date and time with required format and add it to the metadata
    now = datetime.datetime.now()
    date_string = now.strftime("%Y-%m-%d")
    file_header.append("##fileDate={}".format(date_string))
    # Loop through each header lines
    for row in input_file:
        if row.startswith('##'):
            if (row.startswith('##contig') or
                    row.startswith('##FORMAT') or
                    row.startswith('##INFO') or
                    row.startswith('##FILTER') or
                    row.startswith('##sort_order') or
                    row.startswith('##phasing')) and \
                    row.startswith('##INFO=<ID=END') is False:
                file_header.append(row.strip())
        elif row.startswith('#CHROM'):
            file_colnames = row.strip().split('\t')
        else:
            break
    # Get the command-line arguments as a string
    vcf_filter_command = ' '.join(sys.argv)
    # Add BAM_to_gVCF to the metadata
    file_header.append("##source={}".format("VCF_filter/v1.0"))
    file_header.append("##commandline={}".format(vcf_filter_command))
    # write header and metadata
    for line in file_header:
        output_file.write(line + "\n")
    if 'file_colnames' in locals():
        output_file.write("\t".join(file_colnames) + "\n")
    else:
        sys.exit("Column names could not be retrieved")
    # Return
    return file_colnames


# Function to get phasing set sizes
def get_ps_sizes(input_file):
    # Create an empty dictionary
    ps_size_dict = {}
    # Loop over vcf file
    for line in input_file:
        # test if metadata or header
        if not line.startswith('#'):
            # Split into columns and get INFO
            info_col = line.strip().split('\t')[7]
            # Split the info filed into a list based on ";"
            info_col_list = info_col.split(";")
            # Extract only ps information by iterating over the list elements
            ps = ""
            for element in info_col_list:
                # Check if the element starts with "PS="
                if element.startswith("PS="):
                    # Extract the substring starting from "PS="
                    ps = element[3:]
                    break
            # Test if any phasing set information
            if ps == "":
                sys.exit("A PS field is required in INFO for all variants")
            # Count the number occurrences of the PS
            if ps in ps_size_dict:
                ps_size_dict[ps] += 1
            else:
                ps_size_dict[ps] = 1
    return ps_size_dict


# Function to process filtering of a phasing set
def filter_hap(ps_vcf, min_gt_depth, max_gt_depth, min_hap_depth, min_hap_freq):
    # Create two dictionaries
    test_dict = defaultdict(dict)
    ind_dict = defaultdict(dict)
    # Iterate over all individuals to identify and test haplotypes
    for ind in ps_vcf.columns[9:]:
        # Split VCF datapoints
        elt_split = [elt.split(':') for elt in ps_vcf[ind]]
        # Get set of haplotypes by position of the phasing set
        hap_set_by_pos = [elt[0] for elt in elt_split]
        # Test if there is any missing value
        if any("." not in item for item in hap_set_by_pos):
            # Split set of haplotypes by "|"
            hap_by_pos = [elt.split('|') for elt in hap_set_by_pos]
            # Generate a dataframe with haplotypes as rows
            # Get AD information
            ad_set_by_pos = [elt[2] for elt in elt_split]
            ad_by_pos = [elt.split(',') for elt in ad_set_by_pos]
            # Generate a dataframe with AD as rows
            # Transform list of alleles ids into str:
            # [['0', '1'], ['0', '0'], ['0', '0'], ['1', '0']] --> ['0001', '1000'], and store into a list
            hap_list = ["".join(x) for x in zip(*hap_by_pos)]
            # Store AD of each haplotype into a list
            ad_list = [int(elt) for elt in ad_by_pos[0]]
            # Get the value of DP for the phasing set
            dp_ps = int(elt_split[0][1])
            # Fill individual dictionary that stores alleles and depth
            ind_dict[ind] = {"AD": ad_list, "GT": hap_by_pos, "KEY": hap_list}
            # Fill test dictionary by ind using filtering parameters
            test_dict[ind] = {hap_list[i]: f'{False if ad_list[i] < min_hap_depth or dp_ps < min_hap_freq else True}'
                              for i in range(0, len(ad_list))}
    # Aggregate dictionary into a dataframe with haplotypes as rows and individuals as columns,
    # values are boolean stating if the haplotype passed filters for the individual
    hap_test_df = pd.DataFrame(test_dict).fillna(False)
    # Add haplotype to be removed in list if only False for a given haplotype row of hap_test_df
    hap_to_remove = [index for index in hap_test_df.index if "True" not in list(hap_test_df.loc[index])]
    # Update the vcf based on selected alleles and filter for read depth
    ps_vcf.iloc[:, 9:] = ps_vcf.iloc[:, 9:].apply(lambda x: update_hap(x, hap_to_remove, ind_dict,
                                                                       min_gt_depth, max_gt_depth),
                                                  result_type='broadcast')
    # Return
    return ps_vcf


# Function to process filtering of a phasing set
def update_hap(ind_value, hap_to_remove, ind_dict, min_gt_depth, max_gt_depth):
    # Test if individual is in dictionary
    if ind_value.name in ind_dict.keys():
        # Get indices of hapotypes to be removed
        index_to_remove = [i for i in range(0, len(ind_dict[ind_value.name]["KEY"])) if
                           ind_dict[ind_value.name]["KEY"][i] in hap_to_remove]
        # Update dictionary
        ind_dict[ind_value.name]["AD"] = [ind_dict[ind_value.name]["AD"][i] for
                                          i in range(len(ind_dict[ind_value.name]["AD"])) if
                                          ind_dict[ind_value.name]["KEY"][i] not in hap_to_remove]
        ind_dict[ind_value.name]["GT"] = [[elem for i, elem in enumerate(inner_lst) if i not in index_to_remove] for
                                          inner_lst in ind_dict[ind_value.name]["GT"]]
        ind_dict[ind_value.name]["KEY"] = [ind_dict[ind_value.name]["KEY"][i] for
                                           i in range(len(ind_dict[ind_value.name]["KEY"])) if
                                           ind_dict[ind_value.name]["KEY"][i] not in hap_to_remove]
        # Build new haplotypes in GT
        ind_dict[ind_value.name]["GT"] = ['|'.join(elt) for elt in ind_dict[ind_value.name]["GT"]]
        # Compute new DP
        dp_ps_new = 0
        for a in range(len(ind_dict[ind_value.name]["AD"])):
            dp_ps_new += ind_dict[ind_value.name]["AD"][a]
            ind_dict[ind_value.name]["AD"][a] = str(ind_dict[ind_value.name]["AD"][a])
        # Test if DP is within DP bounds and build new GT
        if min_gt_depth <= dp_ps_new <= max_gt_depth:
            ind_dict[ind_value.name]["AD"] = ','.join(ind_dict[ind_value.name]["AD"])
            ind_value = [f'{ind_dict[ind_value.name]["GT"][i]}:{dp_ps_new}:{ind_dict[ind_value.name]["AD"]}' for
                         i in range(0, len(ind_dict[ind_value.name]["GT"]))]
        else:
            ind_value = [f'.:.:.' for _ in range(0, len(ind_dict[ind_value.name]["GT"]))]
    else:
        # Build missing GT
        ind_value = [f'.:.:.' for _ in range(0, len(ind_value))]
    # Return
    return ind_value


# Function to reset ALT alleles according to filtered GT
def filter_reset_variant(row_values):
    # Create ALT and ind dictionary
    allele_dict = defaultdict(lambda: 0)
    ind_dict = defaultdict(dict)
    # Test if the 'ALT' column is '.' or '<NON_REF>'
    if row_values.loc["ALT"] == '.' or row_values.loc["ALT"] == '<NON_REF>':
        return None
    # Iterate over the samples
    for ind in range(9, len(row_values)):
        # Split GT
        gt_split = row_values.iloc[ind].split(':')
        # Get alleles ids
        ind_gt = gt_split[0].split('|')
        # Fill alt dictionary
        for i in ind_gt:
            if i != ".":
                allele_dict[int(i)] += 1
        # Fill ind dictionary
        ind_dict[ind] = {"GT": ind_gt, "DP": gt_split[1], "AD": gt_split[2]}
    # Test if only 0 as allele ids or if only missing values
    if (len(allele_dict) == 1 and 0 in allele_dict) or len(allele_dict) == 0:
        return None
    else:
        # List of sorted alt keys
        allele_dict_keys_sorted = sorted(allele_dict.keys())
        # Create conversion dictionary
        conversion_dict = {allele_dict_keys_sorted[i]: i for i in range(0, len(allele_dict.keys()))}
        # Get ALT alleles
        alt_alleles = row_values.loc['ALT'].split(",")
        # Test if the number of alleles in ALT corresponds to the number of keys in allele dictionary
        if 0 in allele_dict and len(allele_dict) == len(alt_alleles) + 1:
            return row_values
        elif 0 not in allele_dict and len(allele_dict) == len(alt_alleles):
            return row_values
        # Get ALT alleles ids
        alt_ids = [i for i in range(1, len(alt_alleles)+1)]
        alt_alleles_new = alt_alleles.copy()
        # loop over ALT alleles
        for i in alt_ids:
            # Test if key is in dictionary
            if i not in allele_dict:
                alt_alleles_new.remove(alt_alleles[i-1])
        # Append ALT
        row_values.loc['ALT'] = ",".join(alt_alleles_new)
        # Loop over ind
        for ind in range(9, len(row_values)):
            # Append GT
            if "." not in ind_dict[ind]["GT"]:
                gt_new = [str(conversion_dict[int(i)]) for i in ind_dict[ind]["GT"]]
            else:
                gt_new = "."
            row_values.iloc[ind] = f'{"|".join(gt_new)}:{ind_dict[ind]["DP"]}:{ind_dict[ind]["AD"]}'
    return row_values


# Apply program
if __name__ == '__main__':
    main()
