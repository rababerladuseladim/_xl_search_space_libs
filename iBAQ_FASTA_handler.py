import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import re
import os

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class IbaqExtraction:
    """takes a MaxQuant results file and
    returns a sorted list of tuples,
    each tuple containing the protein identifier and the respective iBAQ value of the protein"""
    def __init__(self, filename, keep_contaminants=False, index="Protein IDs"):
        self.filename = filename
        self.keep_contaminants = keep_contaminants
        # pandas data frame that holds protein IDs and their iBAQ intensities, absolute and normalized
        self.results = self.read_file(index)

    def read_file(self, index):
        # read MaxQuant result file with protein ID as index
        allowed_indices = ['Protein IDs', 'Majority protein IDs']
        assert index in allowed_indices, "index needs to be in {} but is {}".format(allowed_indices, index)
        raw_file = pd.read_table(self.filename)
        raw_file.index = raw_file[index]
        raw_file[['Potential contaminant', "Reverse"]] = raw_file[['Potential contaminant', "Reverse"]].astype(str)

        # filter out contaminants
        if not self.keep_contaminants:
            raw_file = raw_file[raw_file["Potential contaminant"] != "+"]

        # filter out decoys
        raw_file_wo_reverse = raw_file[raw_file["Reverse"] != "+"]

        unsorted_results = pd.DataFrame(raw_file_wo_reverse.loc[:, ("iBAQ", "Majority protein IDs")])
        sorted_results = unsorted_results.sort_values(by='iBAQ', ascending=False)
        # add column with logarithmic iBAQ values, natural logarithm
        sorted_results['logIBAQ'] = np.log(sorted_results['iBAQ'])
        # add columns with normalized iBAQ values
        sorted_results['normIBAQ'] = sorted_results['iBAQ'] / sorted_results['iBAQ'].max()
        sorted_results['lognormIBAQ'] = np.log(sorted_results['normIBAQ'])
        return sorted_results

    @staticmethod
    def split_up_protein_groups(list_to_clean):
        cleaned_list_of_proteins = []
        for element in list_to_clean:
            subelements = element.split(';')
            # removal of "CON__"
            # subelements = [subelem.split('CON__')[-1] for subelem in subelements]
            # removal of special ":" chars and the characters preceeding this
            # subelements.extend([subelem.split(':')[-1] for subelem in subelements])
            cleaned_list_of_proteins.extend(subelements)
        return cleaned_list_of_proteins

    def get_rel_log_higher_than(self, fraction):
        """returns a list of proteins of at least 'fraction' log(intensity/maximum iBAQ intensity)"""
        if np.isneginf(fraction):   # to include '-inf' log values
            list_of_proteins = self.results.index.tolist()
        else:
            top_frac = self.results[self.results['lognormIBAQ'] >= fraction]
            list_of_proteins = top_frac.index.tolist()
        # split single elements containing multiple proteins into multiple elements
        cleaned_list_of_proteins = self.split_up_protein_groups(list_of_proteins)
        return cleaned_list_of_proteins

    def get_rel_higher_than(self, fraction):
        """returns a list of proteins of at least 'percentage' intensity of the maximum iBAQ value"""
        if fraction == 0:   # to include '-inf' log values
            list_of_proteins = self.results.index.tolist()
        else:
            top_frac = self.results[self.results['normIBAQ'] >= fraction]
            list_of_proteins = top_frac.index.tolist()
        # split single elements containing multiple proteins into multiple elements
        cleaned_list_of_proteins = self.split_up_protein_groups(list_of_proteins)
        return cleaned_list_of_proteins

    def get_top_quant(self, quantile):
        """returns a list of proteins that have an intensity >= the specified quantile"""
        top_quant = self.results[self.results['iBAQ'] >= self.results.quantile(quantile)['iBAQ']]
        top_quant = top_quant.index.tolist()
        top_quant_list = self.split_up_protein_groups(top_quant)
        return top_quant_list

    def get_top_no(self, number_of_proteins):
        """teakes the first "number_of_proteins" proteins from self.result.index"""
        top_no = self.results.iloc[0:number_of_proteins, :]
        top_no = top_no.index.tolist()
        top_no_list = self.split_up_protein_groups(top_no)[0:number_of_proteins]
        return top_no_list

    def visualization(self):
        """histogram comparison of the different preprocessed iBAQ intensities"""
        bins = 100
        vals = self.results['normIBAQ'].dropna()
        logvals = self.results[self.results['normlogIBAQ'] >= 0]['normlogIBAQ']

        plt.subplot(211)    # the first subplot in the first figure (figure consists of two rows)
        # x-axis ends at the end of the values
        # plt.xlim(0, 1)
        plt.hist(vals, bins=bins)
        # plt.xlabel('normalized intensity')
        plt.ylabel('Frequency')
        plt.title('normalized intensity')
        plt.grid(True)

        # get 2 rows, 1 column and fignum
        plt.subplot(212)    # the second subplot in the first figure (figure consists of two rows)
        # x-axis ends at the end of the values
        # plt.xlim(0, 1)
        # the histogram of the data
        plt.hist(logvals, bins=bins)
        plt.xlabel('normalized log intensity')
        plt.ylabel('Frequency')
        # plt.title('normalized logarithmic values')
        plt.grid(True)

        plt.show()


# # test cases
# ibaq_list = IbaqExtraction("../../../Data/Input/170323_iBAQ_based_opt/Fr14/MaxQuant_results_Fr14_proteinGroups.txt",
#                            keep_contaminants=True, index="Majority protein IDs")
# print len(ibaq_list.get_rel_higher_than(0.01))
# # for elem in ibaq_list.get_top_no(200):
# #     print elem
# # print "No. of proteins for 'Protein IDs': {}".format(len(ibaq_list.get_top_no(200)))
# #
# ibaq_list.results.index = ibaq_list.results["Majority protein IDs"]
# print "No. of proteins for 'Majority protein IDs': {}".format(len(ibaq_list.get_top_no(200)[0:200]))

# protein_list = ibaq_list.get_intensity_top_log_frac(0.9)
# print ibaq_list.get_intensity_top_log_frac(0)
# print ibaq_list.normalized_results
# dfgui.show(ibaq_list.results)
# ibaq_list.visualization()
# with open("../testing/prots.txt", "w+") as f:
#     f.write("\n".join(ibaq_list.get_top_quant(0.0)))
# print ibaq_list.get_top_quant(0.5)


class FastaHandler:
    def __init__(self, fasta_filename, re_id_pattern=r'^>.*\|(.*)\|.*'):
        self.filename = fasta_filename
        self.dict = self.read_fasta(re_id_pattern)

    def read_fasta(self, re_id_pattern):
        """read self.filename fasta file into a dict with protein id as key and its sequence as value"""
        mydict = {}
        dict_key = ""
        list_of_non_unique_ids = []
        protein_id_regex = re.compile(re_id_pattern)
        protein_seq_regex = re.compile(r'^[A-Za-z]\B')
        with open(self.filename) as db:
            for line in db.readlines():
                protein_id_regex_hit = protein_id_regex.match(line)
                protein_seq_regex_hit = protein_seq_regex.match(line)
                if protein_id_regex_hit:
                    dict_key = protein_id_regex_hit.group(1)
                    if dict_key not in mydict.keys():
                        mydict[dict_key] = ''
                        protein_not_in_mydict = True
                    else:
                        list_of_non_unique_ids.append(dict_key)
                        protein_not_in_mydict = False
                elif protein_seq_regex_hit and protein_not_in_mydict:
                    mydict[dict_key] += line.strip()
        if not list_of_non_unique_ids:
            logging.debug("no duplicates in fasta file '{}'".format(self.filename))
        else:
            msg = "Duplicates in fasta file '{}': \n{}".format(self.filename, list_of_non_unique_ids)
            print msg
            logging.warning(msg)
        return mydict

    def build_fasta(self, protein_id_list, filename, max_number=None):
        bool_unfound_protein = False
        no_found_protein = 0
        path, filename_only = os.path.split(filename)
        filename_wo_ext, ext = os.path.splitext(filename_only)
        filename_not_found = os.path.join(path, filename_wo_ext + "_not_found_in_reference" + ext)
        with open(filename_not_found, 'w+') as f_not_found:
            f_not_found.write(r"# all the proteins not found in the reference fasta are listed here."+'\n')
            with open(filename, 'w+') as f:
                for protein_id in protein_id_list:
                    if protein_id in self.dict.keys():
                        f.write('>' + protein_id + '\n')
                        f.write(self.dict[protein_id] + '\n')
                        no_found_protein += 1
                        if max_number:
                            if no_found_protein == max_number:
                                break
                    else:
                        f_not_found.write('>' + protein_id + '\n')
                        bool_unfound_protein = True
                        # msg = "protein ID '{}' not in fasta file '{}'".format(protein_id, self.filename)
                        # print msg

        if bool_unfound_protein:
            logging.warning("some of the specified proteins were not found in {0}, please check {1}"
                            .format(self.filename, filename_not_found))

# # # test cases
if __name__ == "__main__":
    ibaq_list = IbaqExtraction("../../../Data/Input/170323_iBAQ_based_opt/170704-reviewed_fasta/Fr14/proteinGroups.txt",
                               keep_contaminants=True, index="Majority protein IDs")
    # fasta_object = FastaHandler(r"../../../Data/Input/170323_iBAQ_based_opt/fasta_HomoSapiens_reviewed_UP000005640_170630/uniprot-reviewed-yes+AND+proteome-up000005640.fasta")
    # fasta_object.build_fasta(ibaq_list.get_rel_higher_than(0.01), r"../test/fasta.fasta")
