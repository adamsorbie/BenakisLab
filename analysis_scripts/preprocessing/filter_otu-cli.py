import pandas as pd
import argparse
import sys
import os


# to-do
# 1. comment code
# 2. add useful help and examples
# 3. rename arguments and vars to be more interpretable
# 4. test code
# 5. error handling
# 6. print to console

def read_files(metadata, otu_table, filter_file=None):
    """

    :param metadata:
    :param otu_table:
    :param filter_file:
    :return:
    """
    meta = pd.read_csv(metadata, index_col=0, header=0, sep="\t")
    otu = pd.read_csv(otu_table, index_col=0, header=0, sep="\t")
    if filter_file:
        filter_file = pd.read_csv(filter_file, index_col=0, header=0, sep="\t", names=["filter"])
        return meta, otu, filter_file
    return meta, otu


def write_tab(file1, filename1, file2=None, filename2=None):
    """

    :param file1:
    :param filename1:
    :param file2:
    :param filename2:
    :return:
    """
    file1.to_csv(filename1, sep="\t", encoding="utf-8")
    if file2 is not None:
        # need to also make sure filename2 exists here
        file2.to_csv(filename2, sep="\t", encoding="utf-8")


def make_dir(dirnames):
    """
    :param list dirnames: list of directories to create
    :return:
    """
    for i in dirnames:
        try:
            os.mkdir(i)
        except FileExistsError:
            print("Directory ", i, " already exists")


def filter_by_criteria(metadata, otu_table, col1, filter_by, meta_out=None, otu_out=None, col2=None, filter_by2=None):
    """

    :param metadata:
    :param otu_table:
    :param col1:
    :param filter_by:
    :param meta_out:
    :param otu_out:
    :param col2:
    :param filter_by2:
    :return:
    """
    meta, otu = read_files(metadata, otu_table)
    otu = otu.T
    filtered_meta = meta[meta[col1] == filter_by]
    list_samples = filtered_meta.index.tolist()
    filtered_otu = otu[otu.index.isin(list_samples)].T

    if col2:
        filtered_meta = meta[(meta[col1] == filter_by) & (meta[col2] == filter_by2)]
        list_samples = filtered_meta.index.tolist()
        filtered_otu = otu[otu.index.isin(list_samples)].T
        return filtered_meta, filtered_otu
    # should include check if a) type is string and b) string exists in input df
    if all([meta_out, otu_out]):
        write_tab(filtered_meta, meta_out, filtered_otu, otu_out)


def filter_by_list(sample_or_otu, filter_file, metadata, otu_table, meta_out, otu_out, prefilter=False, col1=None,
                filter_by=None, col2=None, filter_by2=None, ):
    """

    :param sample_or_otu:
    :param filter_file:
    :param metadata:
    :param otu_table:
    :param meta_out:
    :param otu_out:
    :param prefilter:
    :param col1:
    :param filter_by:
    :param col2:
    :param filter_by2:
    :return:
    """
    meta, otu, filter_file = read_files(metadata, otu_table, filter_file=filter_file)
    filter_list = filter_file.index.tolist()
    if prefilter:
        meta, otu = filter_by_criteria(metadata, otu_table, col1=col1, filter_by=filter_by)
        if col2:
            meta, otu = filter_by_criteria(metadata, otu_table, col1=col1, filter_by=filter_by, col2=col2,
                                        filter_by2=filter_by2)
    if sample_or_otu == "sample_list":
        filtered_meta = meta[meta.index.isin(filter_list)]
        filtered_otu = otu[filter_list]
        write_tab(filtered_meta, meta_out, filtered_otu, otu_out)
    elif sample_or_otu == "OTU_list":
        filtered_otu = otu[filter_list]
        write_tab(meta, meta_out, filtered_otu, otu_out)


def filter_combined(combined_otu, col1, filter_by, combined_out, col2=None, filter_by2=None):
    """

    :param combined_otu:
    :param col1:
    :param filter_by:
    :param combined_out:
    :param col2:
    :param filter_by2:
    :return:
    """
    combined_filt = combined_otu[combined_otu[col1] == filter_by]
    if col2:
        combined_filt_2 = combined_filt[combined_filt[col2] == filter_by2]
        write_tab(combined_filt_2, combined_out)
    write_tab(combined_filt, combined_out)


def main():
    """
    parse cli arguments and write output to file
    :return:
    """

    parser = argparse.ArgumentParser(description="DESCRIPTION\n"
                                                "Filter an OTU table or combined OTU table by metadata or filter using"
                                                "list of samples/OTUs\n"
                                                "\n\n==========================BASIC USAGE==========================\n"
                                                "\n$ filter_otu-cli.py --filt_vars True -i OTU.tab -m meta.tab -oo "
                                                "otu_filt.tab -mo meta-filt.tab -c1 phenotype -f sick\n"
                                                , formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--filt_vars", type=bool, help="filter by variable")
    parser.add_argument("--filt_list", type=bool, help="filter by list")
    parser.add_argument("--filt_combined", type=bool, help="filter combined OTU/tax table")
    parser.add_argument("--filt_file", type=str, help="filter from inpur file")
    parser.add_argument("-i", "--input", type=str, help="Input otu table")
    parser.add_argument("-m", "--metadata", type=str, help="metadata file")
    parser.add_argument("-oo", "--otu_output", type=str, help="OTU out")
    parser.add_argument("-mo", "--meta_output", type=str, help="metadata out")
    parser.add_argument("-c1", "--column1", type=str, help="1st column to use for filtering")
    parser.add_argument("-c2", "--column2", type=str, help="2nd column to use for filtering")
    parser.add_argument("-f", "--filter_by", type=str, help="variable to filter by c1")
    parser.add_argument("-f2", "--filter_by2", type=str, help="variable to filter by c2")
    parser.add_argument("--list_filt_type", type=str,
                        help="Filter by sample or OTU. Options are 'sample_list' or 'OTU_list'.")
    parser.add_argument("-p", "--prefilt", type=str, help="if filtering by list, prefilter by variables in addition")

    args = parser.parse_args()


    if args.filt_list:
        filter_by_list(args.list_filt_type, args.filt_file, args.metadata, args.input, args.meta_output, args.otu_output,
                    args.prefilt, args.column1, args.filter_by, args.column2, args.filter_by2)

    elif args.filt_combined:
        filter_combined(args.input, args.column1, args.filter_by, args.otu_output, args.column2, args.filter_by2)

    elif args.filt_vars:
        filter_by_criteria(args.metadata, args.input, args.column1, args.filter_by, args.meta_output, args.otu_output,
                            col2=args.column2, filter_by2=args.filter_by2)
    else:
        print("Error no filtering method selected")
        return None


if __name__ == '__main__':
    main()
