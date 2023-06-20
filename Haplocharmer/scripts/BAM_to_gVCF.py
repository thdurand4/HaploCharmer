# BAM_to_gVCF

# Import modules
import pysam
import argparse
import datetime
import sys
import gzip
import os


# Function main
def main():
    # Create the parser object
    parser = argparse.ArgumentParser(description='BAM_to_gVCF.py is a python script to count allele combinations '
                                                 'shared on same sequencing reads (small haplotypes) in a BAM or '
                                                 'SAM alignment file and report results into a gVCF file.')
    # Add arguments to the parser
    parser.add_argument('-bam', '--input_bam', type=str, help='Input BAM or SAM file.', required=True)
    parser.add_argument('-bed', '--input_bed', type=str, help='Input BED file.', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output gVCF file prefix.', required=True)
    parser.add_argument('-q', '--mapping_quality', type=float, help='Mapping quality threshold', default=1)
    parser.add_argument('--compress', action='store_true', help='Compress output using gzip')
    parser.add_argument('--report_all_positions', action='store_true', help='Report all positions in gVCF')
    parser.add_argument('-n', '--non_variant_alt', type=str,
                        help='Recode non variant ALT, initially coded as <NON_REF>',
                        default=None)
    # Parse the arguments
    args = parser.parse_args()
    # Apply checks
    output_name, alignment_file_mode, gvcf_write_mode, gvcf_read_mode, gvcf_suffix = process_checks(args)
    # Generate phasing_set dictionary
    phasing_set = phasing_set_dict(args.input_bed)
    # Open BAM file
    with pysam.AlignmentFile(args.input_bam, alignment_file_mode) as bam_file:
        # Generate header
        header = gvcf_header(bam_file)
        # Extract sample name
        sample = bam_file.header.get("RG")[0]["SM"]
        # Open gVCF file
        with pysam.VariantFile(output_name + gvcf_suffix, "w", header=header) as output_vcf:
            # Loop through each set of variants
            for phasing_set_key, phasing_set_value in phasing_set.items():
                # Count alleles using the count_alleles function
                var_list = count_haplotypes(phasing_set_key, phasing_set_value, bam_file, args.mapping_quality,
                                            args.report_all_positions)
                if type(var_list) == list:
                    for var_elem in var_list:
                        # Create a new gvcf record
                        record = gvcf_record(output_vcf, var_elem, sample)
                        # Add the VariantRecord object to the VariantFile object
                        output_vcf.write(record)
    # Modify non variant ALT and END field
    if args.non_variant_alt is not None or args.report_all_positions:
        modify_gvcf(output_name, args.non_variant_alt, args.compress, args.report_all_positions, gvcf_suffix,
                    gvcf_read_mode, gvcf_write_mode)


# Function to apply checks
def process_checks(args):
    # Test for the presence of files
    if not os.path.isfile(args.input_bam):
        sys.exit("There is no .bam or .sam files at this path")
    if not os.path.isfile(args.input_bed):
        sys.exit("There is no .bed file at this path")
    # Create output directory if it does not exist already
    output_dir = os.path.dirname(args.output)
    if len(output_dir) > 0:
        os.makedirs(output_dir, exist_ok=True)
    # Remove vcf ang gz extensions from output if any
    output_name, output_extension = os.path.splitext(args.output)
    if output_extension != ".gvcf" or output_extension != ".g" or \
            output_extension != ".vcf" or output_extension != ".gz":
        output_name=args.output
    else:
        while output_extension == ".gvcf" or output_extension == ".g" or \
                output_extension == ".vcf" or output_extension == ".gz":
            output_name, output_extension = os.path.splitext(output_name)
    # Test for alignment file type
    if args.input_bam.endswith(".bam"):
        alignment_file_mode = "rb"
    elif args.input_bam.endswith(".sam"):
        alignment_file_mode = "r"
    else:
        sys.exit("Invalid input type, please use .bam or .sam files")
    # Build gVCF suffix
    gvcf_suffix = ".g.vcf"
    # Test if gVCF compression is needed
    if args.compress:
        gvcf_suffix += ".gz"
        gvcf_write_mode = "wt"
        gvcf_read_mode = "rt"
    else:
        gvcf_suffix += ""
        gvcf_write_mode = "w"
        gvcf_read_mode = "r"
    # Return
    return output_name, alignment_file_mode, gvcf_write_mode, gvcf_read_mode, gvcf_suffix


# Function to get unique elements from a list
def unique(list1):
    # initialize a null list
    unique_list = []
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # return list
    return unique_list


# Function to generate a list with chromosome, position, reference allele and depth for each position
def ref_list(read_ref, phasing_set_value, dp):
    # Create an empty list to store the reference allele and the total depth
    variant_info = []
    # Iterate through the reference
    for position in range(phasing_set_value["START"], phasing_set_value["END"] + 1):
        variant_info.append({
            "CHROM": phasing_set_value["CHROM"],
            "POS": position,
            "REF": read_ref.get_reference_sequence()[position - read_ref.reference_start].upper(),
            "DP": dp
        })
    return variant_info


# Function to generate dictionary to store the haplotype counts with allele and position information
def hap_dict(spanning_reads, phasing_set_value):
    # Create an empty dictionary
    allele_counts = {}
    # Iterate through the reads in the BAM file
    for read in spanning_reads:
        # Initialize a key to store the combination of alleles for this read
        key = ''
        # Soft clip bool to indicate if soft clip is over
        soft_clip_out = False
        # Initialize a list to store the allele information for this read
        alleles = []
        # read alignment information
        read_align_info = read.get_aligned_pairs()  # optional: with_seq=True
        # Current position on the reference
        current_ref_pos = read.reference_start
        # Iterate through the read
        for read_pos, ref_pos in read_align_info:
            # Check if there is an insertion
            if ref_pos is None:
                ref_pos = current_ref_pos
            else:
                current_ref_pos += 1
                soft_clip_out = True
            # Check if the position to be evaluated is within the interval of interest
            if phasing_set_value["START"] <= ref_pos <= phasing_set_value["END"] and soft_clip_out:
                # Check if there is a deletion
                if read_pos is None:
                    allele = "*"
                else:
                    allele = read.seq[read_pos]
                key += allele
                alleles.append({'allele': allele, 'pos': ref_pos})
            elif ref_pos > phasing_set_value["END"]:
                # Exit loop over CIGAR operations
                break
        # Check if this combination of alleles has been seen before
        if key in allele_counts:
            # If it has, increment the count for this combination of alleles
            allele_counts[key]['count'] += 1
        else:
            # If it has not, create a new entry for this combination of alleles and initialize the count to 1
            allele_counts[key] = {'count': 1, 'alleles': alleles}
    return allele_counts


# Function to generate an output list for the phasing set to generate gVCF
def output_gvcf_list(allele_counts, variant_info, phasing_set_key, report_all_positions):
    # Create a list
    output_list = []
    # Non variant incrementation variable
    non_var_incr = 0
    # Iterate through variants
    for variant in variant_info:
        # List of alleles
        allele_nucl = []
        # List of allele depth
        allele_depth = []
        # Iterate through haplotype allele combinations
        for haplotype in allele_counts:
            allele_length = 0
            allele_depth.append(allele_counts[haplotype]["count"])
            # Iterate through positions
            for hap_pos in allele_counts[haplotype]["alleles"]:
                if hap_pos["pos"] < variant['POS']:
                    continue
                elif hap_pos["pos"] == variant['POS'] and allele_length == 0:
                    allele_nucl.append(hap_pos["allele"])
                    allele_length += 1
                elif hap_pos["pos"] == variant['POS'] and allele_length != 0:
                    allele_nucl[-1] += hap_pos["allele"]
                else:
                    break
        # Get unique alleles
        hap_variant_total = [variant['REF']]
        for i in allele_nucl:
            hap_variant_total.append(i)
        hap_variant_unique = unique(hap_variant_total)
        # create dictionary with keys as elements of hap_variant_unique and values as their corresponding indices
        d = {v: i for i, v in enumerate(hap_variant_unique)}
        # map each element of hap_variant to its corresponding index using the dictionary
        hap_variant_id = [d[x] for x in allele_nucl]
        # Get alt alleles
        hap_variant_alt = [x for x in hap_variant_unique if x != variant['REF']]
        # Create a variant dictionary
        var_elem = {
            "CHROM": variant['CHROM'],
            "POS": variant['POS'] + 1,
            "REF": variant['REF'],
            "ALT": hap_variant_alt,
            "INFO": {"PS": phasing_set_key},
            "SAMPLE": {"GT": hap_variant_id, "DP": variant['DP'], "AD": allele_depth}
        }
        # Append list if variant or if beginning of block or if all positions need to be reported
        if len(var_elem["ALT"]) > 0 or report_all_positions:
            var_elem["INFO"]["END"] = var_elem["POS"]
            output_list.append(var_elem)
            non_var_incr = 0
        elif len(var_elem["ALT"]) == 0 and non_var_incr == 0:
            var_elem["INFO"]["END"] = var_elem["POS"]
            output_list.append(var_elem)
            non_var_incr += 1
        else:
            output_list[-1]["INFO"]["END"] += 1
    return output_list


# Function to count haplotypes
def count_haplotypes(phasing_set_key, phasing_set_value, bam_file, mapping_quality, report_all_positions):
    # Fetch the reads that overlap the variant positions
    overlapping_reads = bam_file.fetch(phasing_set_value["CHROM"], phasing_set_value["START"], phasing_set_value["END"])
    # Filter out the reads that don't span the entire set of positions
    spanning_reads = [read for read in overlapping_reads if read.reference_start <= phasing_set_value["START"] and
                      read.reference_end > phasing_set_value["END"] and read.mapping_quality >= mapping_quality]
    # Check if there are no spanning reads, if so exit
    dp = len(spanning_reads)
    if dp == 0:
        print("Warning: No spanning reads found for:", phasing_set_key, "(not considered)")
        return -1
    # Iterate over the reads until you find one with an MD tag to retrieve the reference sequence
    for read_00 in spanning_reads:
        if read_00.has_tag("MD"):
            read_0 = read_00
            break
    # Check if there are no reads with MD tag, if so exit as the reference alleles cannot be retrieved
    if 'read_0' not in locals():
        print("Warning: No reads with MD tag found for:", phasing_set_key, "(not considered)")
        return -1
    # Create a list with chromosome, position, reference allele and depth for each position
    phasing_set_info = ref_list(read_0, phasing_set_value, dp)
    # Delete read_0
    del read_0
    # Dictionary to store the haplotype counts with chromosome, allele and position information
    allele_counts = hap_dict(spanning_reads, phasing_set_value)
    # Create list of gVCF rows as output
    output_list = output_gvcf_list(allele_counts, phasing_set_info, phasing_set_key, report_all_positions)
    return output_list


# Function to create the phasing set dictionary
def phasing_set_dict(bed_file_path):
    with open(bed_file_path, 'r') as bed_file:
        # Create an empty dictionary to store the variants
        variants = {}
        # Loop through each line in the bed file
        for line in bed_file:
            # Split the line into its fields
            fields = line.strip().split()
            # Check that the line contains four fields
            if len(fields) != 4:
                sys.exit(f"Invalid line in bed file: {line}")
            # Extract the chromosome, start position, end position and name
            chromosome = fields[0]
            start = int(fields[1]) - 1
            end = int(fields[2]) - 1
            name = fields[3]
            # Create an empty list to store the variants for this line
            line_variants = {"CHROM": chromosome, "START": start, "END": end}
            # Create a new entry for this line
            variants[name] = line_variants
    return variants


# Function to generate gVCF header
def gvcf_header(bam_file):
    # Create a new header object
    header = pysam.VariantHeader()
    # Get the current date and time with required format and add it to the metadata
    now = datetime.datetime.now()
    date_string = now.strftime("%Y-%m-%d")
    header.add_meta("fileDate", date_string)
    # Get the header of the BAM file
    bam_header = bam_file.header
    # Extract and add the name of the reference sequence to the metadata
    if getattr(bam_file, "reference_filename") is not None:
        header.add_meta("reference", bam_file.reference_filename)
    # Extract and add the contig information to the metadata
    contigs = [{"name": sq["SN"], "length": sq["LN"]} for sq in bam_header.get("SQ")]
    for contig in contigs:
        header.add_line("##contig=<ID={},length={}>".format(contig["name"], contig["length"]))
    # Add phasing information to the metadata
    header.add_meta("phasing", "partial")
    # Extract and add the program information to the metadata
    programs = [{"name": pg["ID"], "version": pg.get("VN", ""), "command": pg.get("CL", "")} for pg in
                bam_header.get("PG")]
    for program in programs:
        header.add_meta("source", "/".join([program["name"], program["version"]]))
        header.add_meta("commandline", program["command"])
    # Get the command-line arguments as a string
    bam_to_gvcf_command = ' '.join(sys.argv)
    # Add BAM_to_gVCF to the metadata
    header.add_meta("source", "/".join(["BAM_to_gVCF", "v1.0"]))
    header.add_meta("commandline", bam_to_gvcf_command)
    # Extract and add the sort order to the metadata
    sort_order = bam_header.get("HD").get("SO", "")
    header.add_meta("sort_order", sort_order)
    # Add INFO and FORMAT to the metadata
    header.add_meta(key="FORMAT",
                    items=[('ID', 'GT'), ('Number', '1'), ('Type', 'String'), ('Description', 'Genotype')])
    header.add_meta(key="FORMAT",
                    items=[('ID', 'DP'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Read Depth')])
    header.add_meta(key="FORMAT", items=[('ID', 'AD'), ('Number', '.'), ('Type', 'Integer'), (
        'Description', 'Allelic depths for the reference and alternative alleles in the order listed')])
    header.add_meta(key="INFO",
                    items=[('ID', 'PS'), ('Number', '1'), ('Type', 'String'), ('Description', 'Phase set ID')])
    # Extract and add the sample information to the metadata
    sample = bam_header.get("RG")[0]["SM"]
    header.add_sample(sample)
    return header


# Function to generate gVCF record
def gvcf_record(output_vcf, var_elem, sample):
    # Create a new record using the new_record() method of the VariantFile object
    record = output_vcf.new_record()
    # Set the attributes of the new record
    record.chrom = var_elem["CHROM"]
    record.pos = int(var_elem["POS"])
    record.stop = int(var_elem["INFO"]["END"])
    record.ref = var_elem["REF"]
    if len(var_elem["ALT"]) > 0:
        record.alts = var_elem["ALT"]
    record.samples[sample]['GT'] = var_elem["SAMPLE"]["GT"]
    record.samples[sample].phased = True
    record.samples[sample]['DP'] = var_elem["SAMPLE"]["DP"]
    record.samples[sample]['AD'] = var_elem["SAMPLE"]["AD"]
    record.info.update(var_elem["INFO"])
    return record


# Function to modify non variant ALT and END field
def modify_gvcf(output, non_variant_alt, compress, report_all_positions, gvcf_suffix, gvcf_read_mode, gvcf_write_mode):
    with (gzip.open if compress else open)(output + gvcf_suffix, gvcf_read_mode) as output_file, (
            gzip.open if compress else open)(output + gvcf_suffix + ".temp", gvcf_write_mode) as output_file_temp:
        # Loop over file lines
        for line in output_file:
            # Test if metadata or header
            if line.startswith('#'):
                # Write line
                output_file_temp.write(line)
            else:
                # Test if non-variant ALT needs to be modified
                if non_variant_alt is not None:
                    # Replace ALT
                    line = line.replace("<NON_REF>", non_variant_alt)
                # Test if all position need to be reported, if so END field is not necessary and thus removed
                if report_all_positions:
                    end_pos = line.find("END")
                    punct_pos = line.find(";")
                    if end_pos > 0:
                        line = line[:end_pos] + line[punct_pos + 1:]
                # Write line
                output_file_temp.write(line)
    # Replace the original gVCF with the temp gVCF
    os.replace(output + gvcf_suffix + ".temp", output + gvcf_suffix)


# Apply program
if __name__ == '__main__':
    main()
