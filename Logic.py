from Bio import SeqIO

import Util
import LogicPrep

class Logics:
    def __init__(self, brcd_list):
        self.tmp = ""
        self.brcd_list = brcd_list

    def get_seq_from_FASTQ(self, brcd_list, fastq_dict):
        result_list = []
        for val_arr in brcd_list:
            for key, fasq_list in fastq_dict.items():
                for fastq_seq in fasq_list:
                    # check barcode
                    if val_arr[1] in fastq_seq:
                        tmp_arr = [val_arr[0], val_arr[1]]
                        # check original seq
                        if val_arr[2] in fastq_seq:
                            tmp_arr.append(val_arr[2])
                        else:
                            tmp_arr.append("")
                        # check edited seq
                        if val_arr[3] in fastq_seq:
                            tmp_arr.append(val_arr[3])
                        else:
                            tmp_arr.append("")
                        # append FASTQ seq
                        tmp_arr.append(fastq_seq)

                        result_list.append(tmp_arr)

        return result_list

    def get_multi_seq_from_FASTQ(self, fastq_list, idx_arr, result_list):
        print("get_multi_seq_from_FASTQ starts " + str(idx_arr[0]) + " ~ " + str(idx_arr[1]))
        for val_arr in self.brcd_list:
            for fastq_seq in fastq_list[idx_arr[0]:idx_arr[1]]:
                # check barcode
                if val_arr[1] in fastq_seq:
                    tmp_arr = [val_arr[0], val_arr[1]]
                    # check original seq
                    if val_arr[2] in fastq_seq:
                        tmp_arr.append(val_arr[2])
                    else:
                        tmp_arr.append("")
                    # check edited seq
                    if val_arr[3] in fastq_seq:
                        tmp_arr.append(val_arr[3])
                    else:
                        tmp_arr.append("")
                    # append FASTQ seq
                    tmp_arr.append(fastq_seq)

                    result_list.append(tmp_arr)
        return

    def get_list_multi_p_seq_from_FASTQ(self, fastq_list):
        result_list = []
        print("get_multi_p_seq_from_FASTQ starts ")
        for val_arr in self.brcd_list:
            for fastq_seq in fastq_list:
                # check barcode
                if val_arr[1] in fastq_seq:
                    tmp_arr = [val_arr[0], val_arr[1]]
                    # check original seq
                    if val_arr[2] in fastq_seq:
                        tmp_arr.append(val_arr[2])
                    else:
                        tmp_arr.append("")
                    # check edited seq
                    if val_arr[3] in fastq_seq:
                        tmp_arr.append(val_arr[3])
                    else:
                        tmp_arr.append("")
                    # append FASTQ seq
                    tmp_arr.append(fastq_seq)

                    result_list.append(tmp_arr)
        return result_list

    def get_dict_multi_p_seq_from_FASTQ(self, fastq_list):
        result_dict = {}
        print("get_dict_multi_p_seq_from_FASTQ starts ")
        for val_arr in self.brcd_list:
            tttt_brcd = val_arr[1]
            ori_seq = val_arr[2]
            edit_seq = val_arr[3]
            for fastq_seq in fastq_list:
                # check barcode
                if tttt_brcd in fastq_seq:
                    if tttt_brcd not in result_dict:
                        result_dict.update({tttt_brcd: {"Original sequence": 0, "Edited sequence": 0}})

                    # check original seq
                    if ori_seq in fastq_seq:
                        result_dict[tttt_brcd]["Original sequence"] += 1

                    # check edited seq
                    if edit_seq in fastq_seq:
                        result_dict[tttt_brcd]["Edited sequence"] += 1

        return result_dict

    def merge_pool_list(self, pool_list):
        mege_dict = {}
        for data_dict in pool_list:
            for barcd_key, val_dict in data_dict.items():
                if barcd_key in mege_dict:
                    mege_dict[barcd_key]["Original sequence"] += val_dict["Original sequence"]
                    mege_dict[barcd_key]["Edited sequence"] += val_dict["Edited sequence"]
                else:
                    mege_dict.update({barcd_key: val_dict})

        return mege_dict












