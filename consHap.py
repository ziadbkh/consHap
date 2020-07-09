import gzip
import os
import argparse
from argparse import RawDescriptionHelpFormatter
import datetime
#import textwrap
#import sys
#import timeit
#import io
from shutil import copyfile
import timeit
from subprocess import check_call

version = "v0.1"
prog_header = '''\


                    ****************************************
                    *      consHAP                         *
                    *      Version: ''' + version + '''                   *
                    *      Date: 4 Jan 2020                *
                    *      Author: Ziad Al Bkhetan         *
                    *      Email: Ziad.albkhetan@gmail.com *
                    ****************************************              
         '''

help_messages={'out_file':'Prefix of output file.'}

help_messages['help'] = "Show this help message and exit. "

help_messages['in_shape'] = "Prefix of haplotype file in SHAPEIT format (haps, sample). " + \
                               "See https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#formats ."

help_messages['in_hapiur'] = "Prefix of haplotype file in Eigenstrat format (HAPI-UR default output format). " + \
                                "See https://code.google.com/archive/p/hapi-ur/downloads ."

help_messages['in_vcf'] = "Prefix of haplotype file in vcf format. "

help_messages['c2shap'] = "Convert input file (Eigenstrat format or HAPI-UR format) to SHAPEIT format." \
                         " --hapiur and --out are mandatory."

help_messages['c2vcf'] = "Convert input file (in SHAPEIT format) to VCF format. --shapeit and --out are mandatory."

help_messages['cons'] = "Consensus haplotype construction. --file_list and --out are mandatory." \
                       " --reference_imputation is optional."

help_messages['file_list'] = "A list of odd number of phased haplotype files (SHAPEIT format) to be used in consensous"\
                            " constructions. These files should contain the same individuals in the same order."

help_messages['reference_imputation'] = "For phased missing SNPs, this option uses the reference tool imputation" \
                                        " as an output for these cases. The reference tool is the first tool " \
                                        "listed in --file_list ."

help_messages['gzipped'] = "Input and output files are compressed (gz)."


def convert_hapiur_to_shapeit(input_file_prefix, out_file_prefix, gzip_file):
    snps_path = input_file_prefix + ".phsnp"
    sample_path = input_file_prefix + ".phind"
    haps_path = input_file_prefix + ".phgeno"
    if gzip_file:
        haps_path = haps_path + ".gz"

    if not os.path.isfile(snps_path):
        print("missing input file: " + snps_path)
        return

    if not os.path.isfile(sample_path):
        print("missing input file: " + snps_path)
        return

    if not os.path.isfile(haps_path):
        print("missing input file: " + snps_path)
        return

    sample_out_file_path = out_file_prefix + ".sample"
    out_haps_file = out_file_prefix + ".haps"

    #print(snps_path)
    #print(sample_path)
    #print(haps_path)
    try:
        # create sample file
        sample_file = open(sample_out_file_path, "w+")
        sample_file.write("ID_1 ID_2 missing\n0 0 0\n")
        missing = "0"
        with open(sample_path, "r") as individuals_file:
            for individual in individuals_file:
                info = individual.split()[0].split(":")
                sample_file.write(info[0] + " " + info[1][:-2] + " " + missing + "\n")
                next(individuals_file)
        sample_file.close()

    except IOError:
        print("Can't write output file: " + sample_out_file_path)
        return

    #read hapiur SNPs
    out_file = open(out_haps_file, "w+")

    if gzip_file:
        haps_file = gzip.open(haps_path, "rt")
        #out_file = gzip.open(out_haps_file, "wt+")
    else:
        haps_file = open(haps_path, "r")

    snps_file = open(snps_path, "r")

    for hap in haps_file:
        hap = list(hap.strip())
        snp_info = snps_file.readline().split()
        out_file.write(snp_info[1] + " " +
                       snp_info[0] + " " +
                       snp_info[3] + " " +
                       snp_info[5] + " " +
                       snp_info[4] + " ")
        out_file.write(" ".join(hap) + "\n")
    out_file.close()
    haps_file.close()
    snps_file.close()
    if gzip_file:
        check_call(['gzip', out_haps_file])


def convert_shapeit_to_vcf(input_file_prefix, out_file_prefix, gzip_file):

    out_file_path = out_file_prefix +".vcf"
    sample_path = input_file_prefix + ".sample"
    haps_path =  input_file_prefix + ".haps"
    if gzip_file:
        #out_file_path = out_file_path  + ".gz"
        haps_path = haps_path + ".gz"
    if not os.path.isfile(sample_path):
        print("missing input file: " + sample_path)
        return

    if not os.path.isfile(haps_path):
        print("missing input file: " + haps_path)
        return

    individuals = []
    with open(sample_path, "r") as in_f:
        in_f.readline()
        in_f.readline()
        for line in in_f:
            individuals.append(line.split()[1])
    # print(individuals)

    out_file = open(out_file_path, "w+")
    if gzip_file:
        #out_file = gzip.open(out_file_path, "wt+")
        hap_file = gzip.open(haps_path, "rt")
    else:
        hap_file = open(haps_path, "r")

    out_file.write("##fileformat=VCFv4.1\n")
    out_file.write("##fileDate=" + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
    out_file.write("##source=consHAP \n")
    out_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n')
    out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(individuals) + "\n")

    for line in hap_file:
        vals = line.split()
        chrom = vals[0]
        #chrom = "21"
        hap = [vals[i] + "|" + vals[i + 1] for i in range(5, len(vals), 2)]
        out_file.write(chrom + "\t" + vals[2] + "\t" + vals[1] + "\t" +
                       vals[3] + "\t" + vals[4] + "\t.\tPASS\t.\tGT\t" +
                       "\t".join(hap) + "\n")
    hap_file.close()
    out_file.close()
    if gzip_file:
        check_call(['gzip', out_file_path])


def construct_consensus_old(reference_tool_imputation, out_file_name, tool_list):
    if reference_tool_imputation:
        print("Reference tool imputation is ON")
    else:
        print("Reference tool imputation is OFF")

    print("aggregated tools", tool_list)
    reference_tool = tool_list[0]
    print("Reference Tools: ", reference_tool)

    all_tools_phased_data_file_reader = {}
    for tool_output_file in tool_list:
        print("Reading tool's output: ", tool_output_file)
        all_tools_phased_data_file_reader[tool_output_file] = \
            gzip.open(tool_output_file + ".haps.gz", "rt")

    print("Reading individuals from the reference tool's output: ", tool_list[0])

    out_file = gzip.open(out_file_name + ".haps.gz", "wt+")

    previous_phase = None
    i_line = 0
    for line in all_tools_phased_data_file_reader[reference_tool]:
        all_tools_phased_data = []
        #print(line.split()[:20])
        all_tools_phased_data.append(line.strip().split())
        if previous_phase is None:
            previous_phase = []
            for tool_output_file in tool_list:
                previous_phase.append([1] * ((len(all_tools_phased_data[0]) - 5) // 2))

        for tool_output_file in tool_list[1:]:

            vals = all_tools_phased_data_file_reader[tool_output_file].readline().strip().split()
            if vals[3] == all_tools_phased_data[0][3] and vals[4] == all_tools_phased_data[0][4]:
                all_tools_phased_data.append(vals)
            elif vals[4] == all_tools_phased_data[0][3] and vals[3] == all_tools_phased_data[0][4]:
                vals[3] = all_tools_phased_data[0][3]
                vals[4] = all_tools_phased_data[0][4]
                vals[5:] = ['0' if x == '1' else '1' for x in vals[5:]]
                all_tools_phased_data.append(vals)
            else:
                print (i_line, all_tools_phased_data[0][:5])
                print(i_line, tool_output_file, vals[:5])
                print("something wrong! - different alleles for the same SNP")
                print("-------------------")
        i_line += 1
        #print(len(previous_phase))

        cons_hap = [-1] * len(all_tools_phased_data[0])
        cons_hap[:5] = all_tools_phased_data[0][:5]
        #print(len(all_tools_phased_data[0]))
        #print(cons_hap[:25])
        for snp_itr in range(5, len(all_tools_phased_data[0]), 2):
            allele_dic1 = {}
            allele_dic2 = {}
            phase_count = len(set([min(x[snp_itr], x[snp_itr + 1]) +
                     max(x[snp_itr], x[snp_itr + 1]) for x in all_tools_phased_data]))
            if len(set([x[snp_itr] for x in all_tools_phased_data] +
                       [x[snp_itr + 1] for x in all_tools_phased_data])) == 1: #homozygous SNP: one allele for all
                cons_hap[snp_itr] = all_tools_phased_data[0][snp_itr]
                cons_hap[snp_itr + 1] = all_tools_phased_data[0][snp_itr + 1]

            elif not reference_tool_imputation or phase_count == 1: #heterozygous: two alleles for all tools

                for tool_itr in range(len(all_tools_phased_data)):
                    tool_results = all_tools_phased_data[tool_itr]
                    if previous_phase[tool_itr][(snp_itr - 5) // 2] == 1:
                        tmp_phase = tool_results[snp_itr] + tool_results[snp_itr + 1]
                    else:
                        tmp_phase = tool_results[snp_itr + 1] + tool_results[snp_itr]

                    if tmp_phase[0] in allele_dic1.keys():
                        allele_dic1[tmp_phase[0]] += 1
                    else:
                        allele_dic1[tmp_phase[0]] = 1

                    if tmp_phase[1] in allele_dic2.keys():
                        allele_dic2[tmp_phase[1]] += 1
                    else:
                        allele_dic2[tmp_phase[1]] = 1

                phases = sorted(allele_dic1.items(), key=lambda kv: -kv[1])
                cons_hap[snp_itr] = phases[0][0]

                phases1 = sorted(allele_dic2.items(), key=lambda kv: -kv[1])
                cons_hap[snp_itr + 1] = phases1[0][0]

                if phase_count == 1:
                    for tool_itr in range(len(all_tools_phased_data)):
                        if all_tools_phased_data[tool_itr][snp_itr] == phases[0][0]:
                            previous_phase[tool_itr][(snp_itr - 5) // 2] = 1
                        elif all_tools_phased_data[tool_itr][snp_itr] == phases1[0][0]:
                            previous_phase[tool_itr][(snp_itr - 5) // 2] = 2
            else: #reference_tool_imputation:
                 # use the reference values directly
                 if previous_phase[0][(snp_itr - 5) // 2] == 1:
                     cons_hap[snp_itr] = all_tools_phased_data[0][snp_itr]
                     cons_hap[snp_itr + 1] = all_tools_phased_data[0][snp_itr + 1]
                 else:
                     cons_hap[snp_itr] = all_tools_phased_data[0][snp_itr + 1]
                     cons_hap[snp_itr + 1] = all_tools_phased_data[0][snp_itr]


        out_file.write(" ".join(cons_hap) + "\n")

    out_file.close()
    for tool_output_file in tool_list:
        all_tools_phased_data_file_reader[tool_output_file].close()

    copyfile(tool_list[0] + ".sample", out_file_prefix + ".sample")


def construct_consensus(reference_tool_imputation, gzip_file, out_file_name, tool_list):
    if reference_tool_imputation:
        print("Reference tool imputation is ON")
    else:
        print("Reference tool imputation is OFF")

    print("aggregated tools", tool_list)
    reference_tool = tool_list[0]
    print("Reference Tools: ", reference_tool)

    all_tools_phased_data_file_reader = {}
    for tool_output_file in tool_list:
        print("Reading tool's output: ", tool_output_file)
        if gzip_file:
            all_tools_phased_data_file_reader[tool_output_file] = \
            gzip.open(tool_output_file + ".haps.gz", "rt")
        else:
            all_tools_phased_data_file_reader[tool_output_file] = \
                open(tool_output_file + ".haps", "r")

    print("Reading individuals from the reference tool's output: ", tool_list[0])

    out_file = open(out_file_name + ".haps", "w+")

    previous_phase = None
    i_line = 0
    for line in all_tools_phased_data_file_reader[reference_tool]:
        all_tools_phased_data = []
        #print(line.split()[:20])
        all_tools_phased_data.append(line.strip().split())
        if previous_phase is None:
            previous_phase = []
            for tool_output_file in tool_list:
                previous_phase.append([1] * ((len(all_tools_phased_data[0]) - 5) // 2))

        for tool_output_file in tool_list[1:]:
            vals = all_tools_phased_data_file_reader[tool_output_file].readline().strip().split()
            if vals[3] == all_tools_phased_data[0][3] and vals[4] == all_tools_phased_data[0][4]:
                all_tools_phased_data.append(vals)
            elif vals[4] == all_tools_phased_data[0][3] and vals[3] == all_tools_phased_data[0][4]:
                vals[3] = all_tools_phased_data[0][3]
                vals[4] = all_tools_phased_data[0][4]
                vals[5:] = ['0' if x == '1' else '1' for x in vals[5:]]
                all_tools_phased_data.append(vals)
            else:
                print (i_line, all_tools_phased_data[0][:5])
                print(i_line, tool_output_file, vals[:5])
                print("something wrong! - different alleles for the same SNP")
                print("-------------------")
        i_line += 1
        #print(len(previous_phase))

        cons_hap = [-1] * len(all_tools_phased_data[0])
        cons_hap[:5] = all_tools_phased_data[0][:5]
        #print(len(all_tools_phased_data[0]))
        #print(cons_hap[:25])
        for snp_itr in range(5, len(all_tools_phased_data[0]), 2):
            allele_dic1 = {}
            allele_dic2 = {}
            phase_count = len(set([min(x[snp_itr], x[snp_itr + 1]) +
                     max(x[snp_itr], x[snp_itr + 1]) for x in all_tools_phased_data]))
            if len(set([x[snp_itr] for x in all_tools_phased_data] +
                       [x[snp_itr + 1] for x in all_tools_phased_data])) == 1: #homozygous SNP: one allele for all
                cons_hap[snp_itr] = all_tools_phased_data[0][snp_itr]
                cons_hap[snp_itr + 1] = all_tools_phased_data[0][snp_itr + 1]

            elif not reference_tool_imputation or phase_count == 1: #heterozygous: two alleles for all tools

                for tool_itr in range(len(all_tools_phased_data)):
                    tool_results = all_tools_phased_data[tool_itr]
                    if previous_phase[tool_itr][(snp_itr - 5) // 2] == 1:
                        tmp_phase = tool_results[snp_itr] + tool_results[snp_itr + 1]
                    else:
                        tmp_phase = tool_results[snp_itr + 1] + tool_results[snp_itr]

                    if tmp_phase[0] in allele_dic1.keys():
                        allele_dic1[tmp_phase[0]] += 1
                    else:
                        allele_dic1[tmp_phase[0]] = 1

                    if tmp_phase[1] in allele_dic2.keys():
                        allele_dic2[tmp_phase[1]] += 1
                    else:
                        allele_dic2[tmp_phase[1]] = 1

                phases = sorted(allele_dic1.items(), key=lambda kv: -kv[1])
                cons_hap[snp_itr] = phases[0][0]

                phases1 = sorted(allele_dic2.items(), key=lambda kv: -kv[1])
                cons_hap[snp_itr + 1] = phases1[0][0]

                if phase_count == 1:
                    for tool_itr in range(len(all_tools_phased_data)):
                        if all_tools_phased_data[tool_itr][snp_itr] == phases[0][0]:
                            previous_phase[tool_itr][(snp_itr - 5) // 2] = 1
                        elif all_tools_phased_data[tool_itr][snp_itr] == phases1[0][0]:
                            previous_phase[tool_itr][(snp_itr - 5) // 2] = 2
            else: #reference_tool_imputation:
                 # use the reference values directly
                 if previous_phase[0][(snp_itr - 5) // 2] == 1:
                     cons_hap[snp_itr] = all_tools_phased_data[0][snp_itr]
                     cons_hap[snp_itr + 1] = all_tools_phased_data[0][snp_itr + 1]
                 else:
                     cons_hap[snp_itr] = all_tools_phased_data[0][snp_itr + 1]
                     cons_hap[snp_itr + 1] = all_tools_phased_data[0][snp_itr]


        out_file.write(" ".join(cons_hap) + "\n")

    out_file.close()
    if gzip_file:
        check_call(['gzip', out_file_name + ".haps"])

    for tool_output_file in tool_list:
        all_tools_phased_data_file_reader[tool_output_file].close()

    copyfile(tool_list[0] + ".sample", out_file_prefix + ".sample")


def initialise_argument_parser():

    parser = argparse.ArgumentParser(prog="consHAP",
                                     formatter_class=RawDescriptionHelpFormatter,
                                     #description=textwrap.dedent(prog_header),
                                     add_help=True)

    parser.add_argument("-v", '--version', action='version', version='%(prog)s ' + version)

    parser.add_argument("-O", "--out", dest="out_file", help=help_messages['out_file'])

    parser.add_argument("-SF", "--shapeit", dest="shapeit_file", help=help_messages['in_shape'])

    parser.add_argument("-HF", "--hapiur", dest="hapiur_file", help=help_messages['in_hapiur'])

    parser.add_argument("-GZ", "--gzipped", dest="is_gzipped", action="store_true",
                        default=False, help=help_messages['gzipped'])

    parser.add_argument("-C2S", "--convert2SHAPEIT", action="store_true",
                        default=False, dest="convert2shapeit", help=help_messages['c2shap'])
    parser.add_argument("-C2V", "--convert2VCF", action="store_true",
                        default=False, dest="convert2vcf", help=help_messages['c2vcf'])
    parser.add_argument("-CONS", "--consensus", action="store_true",
                        default=False, dest="consensus", help=help_messages['cons'])

    parser.add_argument("-F", "--file_list", dest="file_list", nargs='+', help=help_messages['file_list'])

    parser.add_argument("-RI", "--reference_imputation", dest="reference_imputation", action="store_true",
                        default=False, help=help_messages['reference_imputation'])

    return parser


def check_hapiur_file(input_file_prefix, gzip_file):
    all_good = True
    snps_path = input_file_prefix + ".phsnp"
    sample_path = input_file_prefix + ".phind"
    haps_path = input_file_prefix + ".phgeno"
    if gzip_file:
        haps_path = haps_path + ".gz"

    if not os.path.isfile(snps_path):
        print(snps_path + ": File not exist!")
        all_good = False

    if not os.path.isfile(sample_path):
        print(sample_path + ": File not exist!")
        all_good = False

    if not os.path.isfile(haps_path):
        print(haps_path + ": File not exist!")
        all_good = False

    return all_good


def check_shapeit_file(input_file_prefix, gzip_file):
    all_good = True
    sample_path = input_file_prefix + ".sample"
    haps_path = input_file_prefix + ".haps"
    if gzip_file:
        haps_path = haps_path + ".gz"

    if not os.path.isfile(sample_path):
        print(sample_path + ": File not exist!")
        all_good = False

    if not os.path.isfile(haps_path):
        print(haps_path + ": File not exist!")
        all_good = False

    return all_good


def check_output_file(file_name):
    if not os.path.isfile(file_name):
        return True
    else:
        print(file_name+ ": File exist!")
        return False


def check_file_list(file_list, gzip_file):
    if len(file_list) % 2 != 1:
        print ("An odd number of files is required for majority voting!")
        return False

    all_good = True
    for f_nm in file_list:
        if not check_shapeit_file(f_nm, gzip_file):
            all_good = False
    return all_good


def print_help(optional_message=""):

    print(optional_message)
    print("Use the command --help or -h to see all available options and parameters!")
    #parser.print_help(sys.stderr)


if __name__ == '__main__':

    print(prog_header)

    parser = initialise_argument_parser()
    args = parser.parse_args()

    commands=0

    if args.convert2shapeit:
        commands += 1

    if args.convert2vcf:
        commands += 1

    if args.consensus:
        commands += 1

    if commands != 1:
        optional_message = "Please select only one command to be performed"
        print_help(optional_message)
        exit()

    # check all mandatory fields
    if not args.out_file:
        optional_message = "output file prefix is missing"
        print_help(optional_message)
        exit()
    else:
        out_file_prefix = args.out_file

    if args.convert2vcf:
        if not args.shapeit_file:
            optional_message = "--shap argument is missing! "
            print_help(optional_message)
            exit()
        else:
            if check_shapeit_file(args.shapeit_file, args.is_gzipped):
                start = timeit.default_timer()
                convert_shapeit_to_vcf(args.shapeit_file, out_file_prefix, args.is_gzipped)
                stop = timeit.default_timer()
                print('Execution time: ' + str(int(stop - start)) + " seconds.")
            else:
                optional_message = "An issue with the provided SHAPEIT phased file : " + args.shapeit_file
                print_help(optional_message)
                exit()

    elif args.convert2shapeit:
        if not args.hapiur_file:
            optional_message = "--hapur argument is missing! "
            print_help(optional_message)
            exit()
        else:
            if check_hapiur_file(args.hapiur_file, args.is_gzipped):
                start = timeit.default_timer()
                convert_hapiur_to_shapeit(args.hapiur_file, out_file_prefix, args.is_gzipped)
                stop = timeit.default_timer()
                print('Execution time: ' + str(int(stop - start)) + " seconds.")

            else:
                optional_message = "An issue with the provided HAPI-UR phased file : " + args.hapiur_file
                print_help(optional_message)
                exit()

    elif args.consensus:
        if not args.file_list:
            optional_message = "--file_list argument is missing! "
            print_help(optional_message)
            exit()
        else:
            if check_file_list(args.file_list, args.is_gzipped):
                start = timeit.default_timer()
                construct_consensus(args.reference_imputation, args.is_gzipped, out_file_prefix, args.file_list)
                stop = timeit.default_timer()
                print('Execution time: ' + str(int(stop - start)) + " seconds.")

            else:
                print_help()
                exit()