# VCF_to HapPresAbs

# Import modules
import pandas as pd
import sys
import os
import gzip
import datetime
from pyfaidx import Fasta
import argparse
from collections import defaultdict


def main():
    # Create the parser object
    parser = argparse.ArgumentParser(description='VCF_to_HapPresAbs.py is a python script to generate a format '
                                                 'indicating haplotype presence/absence (HPA) from a VCF with'
                                                 'phase set information filtered using VCF_filter.py')
    # Add arguments to the parser
    parser.add_argument('-i', '--input', type=str, help='Filtered input VCF file', required=True)
    parser.add_argument('-r', '--reference', type=str, help='Reference fasta genome to extract haplotype sequences',
                        required=True)
    parser.add_argument('-o', '--output', type=str, help='Output HPA file prefix', required=True)
    parser.add_argument('-v', '--hap_var_info', type=str, help='Output haplotype variant information file prefix',
                        required=True)
    parser.add_argument('--compress', action='store_true', help='Compress output using gzip')
    parser.add_argument('--max_abs_freq', type=float, default=0.01,
                        help='Maximum haplotype frequency (AD/DP) to declare a haplotype absent')
    parser.add_argument('--min_pres_depth', type=int, default=3,
                        help='Minimum haplotype depth (AD) to declare an haplotype present conditional on '
                             'min_pres_freq also being satisfied')
    parser.add_argument('--min_pres_freq', type=float, default=0.04,
                        help='Minimum haplotype frequency (AD/DP) to declare an haplotype present conditional on '
                             'min_pres_depth being also satisfied')
    # Parse the arguments
    args = parser.parse_args()
    # Apply checks
    output_name, output_suffix, hap_info_name, hap_info_suffix, input_read_mode, input_read_pd_mode, \
        output_write_mode, output_write_pd_mode = process_checks(args)
    # Open VCF to get the header, get column names and phasing set sizes
    with (gzip.open if args.input.endswith(".vcf.gz") else open)(args.input, input_read_mode) as input_file, (
            gzip.open if args.compress else open)(output_name + output_suffix, output_write_mode) as output_file, (
            gzip.open if args.compress else open)(hap_info_name + hap_info_suffix, output_write_mode) as \
            hap_var_info_file:
        # Generate header and column names
        file_colnames = generate_header(input_file, output_file, hap_var_info_file)
        # Create a dictionary with the dimension of each phasing set
        ps_sizes = get_ps_sizes(input_file)
    # Generate the filtered HapPresAbs format phasing set by phasing set
    generate_hpa_by_ps(args, ps_sizes, file_colnames, output_name, output_suffix, hap_info_name, hap_info_suffix,
                       input_read_pd_mode, output_write_pd_mode)


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
    # Remove .vcf, .txt, .hpa and .gz extensions from output if any
    output_name, output_extension = os.path.splitext(args.output)
    while output_extension == ".vcf" or output_extension == ".hpa" or output_extension == ".txt" or \
            output_extension == ".gz":
        output_name, output_extension = os.path.splitext(output_name)
    # Build VCF suffix
    output_suffix = ".hpa"
    # Remove .vcf, .txt, .hpa and .gz extensions from hap_info if any
    hap_info_name, hap_info_extension = os.path.splitext(args.hap_var_info)
    while hap_info_extension == ".vcf" or hap_info_extension == ".hpa" or hap_info_extension == ".txt" or \
            hap_info_extension == ".gz":
        hap_info_name, hap_info_extension = os.path.splitext(hap_info_name)
    # Build VCF suffix
    hap_info_suffix = ".txt"
    # Test if gVCF compression is needed
    if args.compress:
        output_suffix += ".gz"
        hap_info_suffix += ".gz"
        output_write_mode = "wt"
        output_write_pd_mode = "gzip"
    else:
        output_write_mode = "w"
        output_write_pd_mode = None
    # Return
    return output_name, output_suffix, hap_info_name, hap_info_suffix, input_read_mode, input_read_pd_mode, \
        output_write_mode, output_write_pd_mode


# Function to generate the header
def generate_header(input_file, output_file, hap_var_info_file):
    # Create header list
    file_header = ['##fileformat=HapPresAbs']
    # Get the current date and time with required format and add it to the metadata
    now = datetime.datetime.now()
    date_string = now.strftime("%Y-%m-%d")
    file_header.append("##fileDate={}".format(date_string))
    # Loop through each header lines
    for row in input_file:
        if row.startswith('##'):
            if row.startswith('##contig') or row.startswith('##sort_order') or row.startswith('##phasing'):
                file_header.append(row.strip())
        elif row.startswith('#CHROM'):
            file_colnames = ['#HAPLOTYPE', 'CHROM', "START", 'END', "SEQUENCE", "FORMAT"] + row.strip().split('\t')[9:]
        else:
            break
    # Get the command-line arguments as a string
    command = ' '.join(sys.argv)
    # Add BAM_to_gVCF to the metadata
    file_header.append("##source={}".format("VCF_to_HapPresAbs/v1.0"))
    file_header.append("##commandline={}".format(command))
    # Write header and metadata of HapPresAbs-info
    for l, line in enumerate(file_header):
        if l == 0:
            line += "-info"
        hap_var_info_file.write(line + "\n")
    # Add FORMAT information
    file_header.append('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype: presence (1) or absence (0) '
                       'of the haplotype">')
    file_header.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth for the phase set">')
    file_header.append('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Haplotype/allele depth for the phase set">')
    # Write header and metadata of HapPresAbs
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


# Function to generate the filtered HapPresAbs format phasing set by phasing set
def generate_hpa_by_ps(args, ps_sizes, hpa_colnames, output_name, output_suffix, hap_info_name, hap_info_suffix,
                       input_read_pd_mode, output_write_pd_mode):
    # Open VCF file
    vcf_reader = pd.read_table(args.input, sep="\t", iterator=True, header=None, comment='#',
                               compression=input_read_pd_mode)
    # Get reference sequence
    reference = Fasta(args.reference, mutable=False)
    # Iterate through phasing sets
    ps_incr = 1
    for ps_key, ps_value in ps_sizes.items():
        # Get chunk of VCF corresponding to the phasing set
        ps_vcf = (vcf_reader.get_chunk(ps_value + ps_incr))
        ps_incr = 0
        # Set column names
        ps_vcf.columns = ['#CHROM', 'POS', "ID", 'REF', "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + hpa_colnames[6:]
        # Convert to hpa format and filter
        ps_hpa, ps_info = convert_to_hpa(ps_vcf, reference, hpa_colnames, args.max_abs_freq, args.min_pres_depth,
                                         args.min_pres_freq)
        # Write filtered vcf block
        ps_hpa.to_csv(output_name + output_suffix, sep="\t", index=False, mode="a", header=False,
                      compression=output_write_pd_mode)
        # Write haplotype information for the PS
        ps_info_sep = pd.DataFrame({'COL': [f"##PS={ps_key}"]})
        ps_info_sep.to_csv(hap_info_name + hap_info_suffix, sep="\t", index=False, mode="a", header=False,
                           compression=output_write_pd_mode)
        ps_info.to_csv(hap_info_name + hap_info_suffix, sep="\t", index=False, mode="a", header=True,
                       compression=output_write_pd_mode)


# Function to convert to hpa format and filter
def convert_to_hpa(ps_vcf, reference, hpa_colnames, max_abs_freq, min_pres_depth, min_pres_freq):
    # Split into columns and get INFO
    info_col = ps_vcf['INFO'].iloc[0]
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
    # Get chromosome, start and end positions
    chrom = ps_vcf['#CHROM'].iloc[0]
    start = ps.split("_")[1]
    end = ps.split("_")[2]
    # Format information
    format = "GT:DP:AD"
    # Create a dictionary
    ind_dict = defaultdict(dict)
    # Iterate over all individuals
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
            hap_list = [x for x in zip(*hap_by_pos)]
            # Store AD of each haplotype into a list
            ad_list = [int(elt) for elt in ad_by_pos[0]]
            # Get the value of DP for the phasing set
            dp_ps = int(elt_split[0][1])
            # Fill individual dictionary that stores alleles and depth
            ind_dict[ind] = {"DP": dp_ps, "AD": ad_list, "KEY": ["-".join(x) for x in hap_list], "MISSING": False}
        else:
            ind_dict[ind] = {"MISSING": True}
    # Unique set of haplotypes ids
    unique_hap = list(set().union(*[set(v['KEY']) for k, v in ind_dict.items() if 'KEY' in v]))
    # Create hpa dataframe and fill information
    ps_hpa = pd.DataFrame(columns=hpa_colnames, index=range(len(unique_hap)))
    ps_hpa['#HAPLOTYPE'] = unique_hap
    ps_hpa['CHROM'] = chrom
    ps_hpa['START'] = start
    ps_hpa['END'] = end
    ps_hpa['FORMAT'] = format
    # Build reference haplotype sequence
    hap_ref_seq = reference[chrom][int(start) - 1:int(end)].seq
    # Generate list of ALT alleles and POS
    alt_list = [list(ps_vcf.loc[i, "REF"]) + ps_vcf.loc[i, "ALT"].split(',') for i in ps_vcf.index]
    pos_list = list([ps_vcf.loc[i, "POS"] for i in ps_vcf.index])
    # Create dictionary to replace REF alleles by ALT alleles
    alt_dict = {pos: [a.replace('*', '') for a in alt] for pos, alt in zip(pos_list, alt_list)}
    # Generate list of haplotype sequences
    hap_seq = []
    # loop through haplotypes
    for hap in unique_hap:
        new_seq = list(hap_ref_seq)
        # loop through positions
        for i, pos in enumerate(alt_dict):
            new_seq[int(pos) - int(start)] = alt_dict[pos][int(hap.split("-")[i])]
        new_seq = ''.join(new_seq)
        hap_seq.append(new_seq)
    # Fill SEQUENCE information
    ps_hpa['SEQUENCE'] = hap_seq
    # Fill hpa dataframe
    ps_hpa_filled = ps_hpa.apply(lambda hpa_row: fill_hpa(hpa_row, ind_dict, ps, max_abs_freq,
                                                          min_pres_depth, min_pres_freq),
                                 result_type='broadcast', axis=1)
    # Create a hap info dataframe and fill variant information
    ps_info = pd.DataFrame(columns=["#PS", "CHROM", "POS", "REF", "ALT"], index=range(ps_vcf.shape[0]))
    ps_info['#PS'] = ps
    ps_info['CHROM'] = list(ps_vcf['#CHROM'])
    ps_info['POS'] = list(ps_vcf['POS'])
    ps_info['REF'] = list(ps_vcf['REF'])
    ps_info['ALT'] = list(ps_vcf['ALT'])
    # Loop through haplotypes of the phase set
    for i, hap in enumerate(ps_hpa_filled['#HAPLOTYPE']):
        # hap_split = hap.split("_")
        # ps_info[hap_split[len(hap_split)-1]] = unique_hap[i].split("-")
        ps_info[hap] = unique_hap[i].split("-")
    # Return
    return ps_hpa_filled, ps_info


# Function to fill hpa row
def fill_hpa(hpa_row, ind_dict, ps, max_abs_freq, min_pres_depth, min_pres_freq):
    # Loop over individuals
    for ind in hpa_row.index[6:]:
        # Check if missing value
        if ind_dict[ind]["MISSING"]:
            pres_abs = "."
            dp = "."
            ad = "."
        # Check if haplotype in KEY
        elif hpa_row["#HAPLOTYPE"] in ind_dict[ind]["KEY"]:
            ad = [ind_dict[ind]["AD"][i] for i, hap in enumerate(ind_dict[ind]["KEY"]) if hap == hpa_row["#HAPLOTYPE"]][
                0]
            dp = ind_dict[ind]["DP"]
            # Check presence absence state according to filters
            if int(ad) >= min_pres_depth and int(ad) / int(dp) > min_pres_freq:
                pres_abs = "1"
            elif int(ad) / int(dp) < max_abs_freq:
                pres_abs = "0:"
            else:
                pres_abs = "."
        else:
            pres_abs = 0
            dp = ind_dict[ind]["DP"]
            ad = 0
        # Fill datapoint
        hpa_row[ind] = f'{pres_abs}:{dp}:{ad}'
    # Make final haplotype name
    hpa_row["#HAPLOTYPE"] = f"{ps}_hap{hpa_row.name + 1}"
    return hpa_row


# Apply program
if __name__ == '__main__':
    main()
