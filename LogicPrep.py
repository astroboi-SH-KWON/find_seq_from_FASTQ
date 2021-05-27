from Bio import SeqIO
import re
import math

import Logic
import Util

class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def get_unit_len_n_remainer(self, trgt_list, div_num):
        total_len = len(trgt_list)
        print("total_len : " + str(total_len))
        flt_unit_len = total_len / div_num
        print("flt_unit_len : " + str(flt_unit_len))
        int_unit_len = int(flt_unit_len)
        print("util_len : " + str(int_unit_len))

        rest_len = total_len
        rest_flot = flt_unit_len - int_unit_len
        print("rest_flot : " + str(rest_flot))
        print("")

        for i in range(div_num):
            rest_len -= int_unit_len
            print("rest_len : " + str(rest_len))
        remainer = math.ceil(int((rest_flot * div_num) * 1000) / 1000)
        print("last : " + str(remainer))

        return int_unit_len, remainer

    def make_list_to_dict(self, total_list):
        result_dict = {}
        for val_arr in total_list:
            barcd_key = val_arr[1]
            ori_seq_cnt = int(val_arr[2])
            edited_seq_cnt = int(val_arr[3])
            barcd_key_cnt = int(val_arr[4])

            if barcd_key in result_dict:
                result_dict[barcd_key]["Original sequence"] += ori_seq_cnt
                result_dict[barcd_key]["Edited sequence"] += edited_seq_cnt
                result_dict[barcd_key]["TTTT_Barcode_cnt"] += barcd_key_cnt
            else:
                result_dict.update({barcd_key: {'Original sequence': ori_seq_cnt, 'Edited sequence': edited_seq_cnt, 'TTTT_Barcode_cnt': barcd_key_cnt}})

        return result_dict

    def make_4seq_list_to_dict(self, total_list):
        result_dict = {}
        for val_arr in total_list:
            barcd_key = val_arr[1]
            barcd_key_cnt = int(val_arr[2])
            wout_edit_seq_cnt = int(val_arr[3])
            edited_seq_cnt = int(val_arr[4])
            pos_1_seq_cnt = int(val_arr[5])
            pos_2_seq_cnt = int(val_arr[6])

            if barcd_key in result_dict:
                result_dict[barcd_key]["TTTT_Barcode_cnt"] += barcd_key_cnt
                result_dict[barcd_key]["Target_sequences_without_edit"] += wout_edit_seq_cnt
                result_dict[barcd_key]["Target_sequences_with_edit_complete"] += edited_seq_cnt
                result_dict[barcd_key]["Position_1_only"] += pos_1_seq_cnt
                result_dict[barcd_key]["Position_2_only"] += pos_2_seq_cnt
            else:
                result_dict.update({barcd_key: {"TTTT_Barcode_cnt": barcd_key_cnt,
                                                "Target_sequences_without_edit": wout_edit_seq_cnt,
                                                "Target_sequences_with_edit_complete": edited_seq_cnt,
                                                "Position_1_only": pos_1_seq_cnt, "Position_2_only": pos_2_seq_cnt}})

        return result_dict

    def add_missing_brcd_to_dict(self, brcd_list, result_dict):
        print("st ::: add_missing_brcd_to_dict")
        for val_arr in brcd_list:
            tttt_brcd = val_arr[1].upper()
            if tttt_brcd not in result_dict:
                result_dict.update(
                    {tttt_brcd: {"Original sequence": 0, "Edited sequence": 0, "TTTT_Barcode_cnt": 0}})
        print("DONE ::: add_missing_brcd_to_dict")

    def list_to_dict_by_ele_key(self, brcd_list, idx):
        print('st list_to_dict_by_ele_key')
        result_dict = {}
        for brcd_arr in brcd_list:
            key = brcd_arr[idx].upper()
            if key == '':
                continue

            if key in result_dict:
                print(key, ':::::::: duplicated')
                result_dict[key].append(brcd_arr)
            else:
                result_dict.update({key: [brcd_arr]})

        print('DONE list_to_dict_by_ele_key')
        return result_dict

    def check_brcd_length(self, brcd_dict, brcd_path):
        print("st check_brcd_length")
        cnt = 0
        ini_brcd_len = 0
        for brcd_key in brcd_dict.keys():
            if brcd_key == '':
                print("there is empty key value")
                continue
                # raise Exception
            if cnt == 0:
                ini_brcd_len = len(brcd_key)
                new_brcd_len = len(brcd_key)
            else:
                new_brcd_len = len(brcd_key)

            if ini_brcd_len != new_brcd_len:
                print("check barcode length in input file :", brcd_path)
                exit()
            cnt += 1
        print("DONE check_brcd_length")
        return ini_brcd_len
