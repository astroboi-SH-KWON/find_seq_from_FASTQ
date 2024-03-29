from Bio import SeqIO

import Util
import LogicPrep

class Logics:
    def __init__(self, brcd_obj):
        self.tmp = ""
        self.brcd_obj = brcd_obj
        self.brcd_pos_arr = [0, 0]
        self.brcd_len = 20
        self.rev_com_flag = False

    def complement_char(self, ch):
        complement_char_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        try:
            return complement_char_dict[ch]
        except:
            print("complement_char : [" + ch + "]")
            raise Exception

    def make_complement_string(self, trgt_seq):
        comp_seq = ""
        for ch in trgt_seq:
            try:
                comp_seq += self.complement_char(ch)
            except:
                raise Exception
        return comp_seq

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
        # self.brcd_obj is list
        for val_arr in self.brcd_obj:
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
        print("get_list_multi_p_seq_from_FASTQ starts ")
        # self.brcd_obj is list
        for val_arr in self.brcd_obj:
            index_name = val_arr[0]
            tttt_brcd = val_arr[1].upper()
            ori_seq = val_arr[2].upper()
            edit_seq = val_arr[3].upper()
            for fastq_str in fastq_list:
                fastq_seq = fastq_str.upper()
                # check barcode
                if tttt_brcd in fastq_seq:
                    tmp_arr = [index_name, tttt_brcd]

                    # TODO check seq after tttt_brcd seq
                    fastq_seq_aft_barcd = fastq_seq[fastq_seq.index(tttt_brcd) + len(tttt_brcd):]
                    # check original seq
                    if ori_seq in fastq_seq_aft_barcd:
                        tmp_arr.append(ori_seq)
                    else:
                        tmp_arr.append("")
                    # check edited seq
                    if edit_seq in fastq_seq_aft_barcd:
                        tmp_arr.append(edit_seq)
                    else:
                        tmp_arr.append("")
                    # append FASTQ seq
                    tmp_arr.append(fastq_seq)

                    result_list.append(tmp_arr)
        return result_list

    def check_barcode_freq(self, brcd_arr, fastq_seq, result_dict):
        tttt_brcd = brcd_arr[1].upper()
        ori_seq = brcd_arr[2].upper()
        edit_seq = brcd_arr[3].upper()

        if tttt_brcd not in result_dict:
            result_dict.update(
                {tttt_brcd: {"Original sequence": 0, "Edited sequence": 0, "TTTT_Barcode_cnt": 1}})
        else:
            result_dict[tttt_brcd]["TTTT_Barcode_cnt"] += 1

        # check seq after tttt_brcd seq
        fastq_seq_aft_barcd = fastq_seq[fastq_seq.index(tttt_brcd) + len(tttt_brcd):]
        # check original seq
        if ori_seq in fastq_seq_aft_barcd:
            result_dict[tttt_brcd]["Original sequence"] += 1

        # check edited seq
        if edit_seq in fastq_seq_aft_barcd:
            result_dict[tttt_brcd]["Edited sequence"] += 1

    def get_dict_by_brcd_w_brcd_range_from_FASTQ(self, fastq_list):
        result_dict = {}
        print("st get_dict_by_brcd_w_brcd_range_from_FASTQ")
        # self.brcd_obj is dict
        brcd_dict = self.brcd_obj

        for fastq_str in fastq_list:
            fastq_seq = fastq_str.upper()
            if self.rev_com_flag:
                fastq_seq = self.make_complement_string(fastq_str.upper())[::-1]

            en_brcd_i = self.brcd_pos_arr[1]
            if en_brcd_i < 0:
                en_brcd_i += len(fastq_seq)
                if en_brcd_i + self.brcd_len > len(fastq_seq):
                    en_brcd_i = len(fastq_seq) - self.brcd_len
            elif en_brcd_i == 0 or en_brcd_i + self.brcd_len > len(fastq_seq):
                en_brcd_i = len(fastq_seq) - self.brcd_len

            for brcd_i in range(self.brcd_pos_arr[0], en_brcd_i):
                candi_brcd = fastq_seq[brcd_i: brcd_i + self.brcd_len]

                if candi_brcd in brcd_dict:
                    self.check_barcode_freq(brcd_dict[candi_brcd][0], fastq_seq, result_dict)

        print("DONE get_dict_by_brcd_w_brcd_range_from_FASTQ")
        return result_dict

    def check_barcode_freq_w_nSeq(self, brcd_arr, fastq_seq, result_dict):

        tttt_brcd = ""
        fastq_seq_aft_barcd = ""
        for i in range(len(brcd_arr)):
            tmp_seq = brcd_arr[i].upper()
            # # skip INDEX
            if i == 0 or tmp_seq.replace(" ", "") == '':
                continue

            # # check brcd
            if i == 1:
                tttt_brcd = tmp_seq
                if tttt_brcd not in result_dict:
                    result_dict.update({tttt_brcd: {"Original sequence": 0, "TTTT_Barcode_cnt": 1}})
                else:
                    result_dict[tttt_brcd]["TTTT_Barcode_cnt"] += 1
                # check seq after tttt_brcd seq
                fastq_seq_aft_barcd = fastq_seq[fastq_seq.index(tttt_brcd) + len(tttt_brcd):]

            # # check original seq
            elif i == 2:
                ori_seq = tmp_seq
                if ori_seq in fastq_seq_aft_barcd:
                    result_dict[tttt_brcd]["Original sequence"] += 1

            # # check edited seq
            else:
                edit_seq = tmp_seq
                if edit_seq in fastq_seq_aft_barcd:
                    if str(i) in result_dict[tttt_brcd]:
                        result_dict[tttt_brcd][str(i)] += 1
                    else:
                        result_dict[tttt_brcd].update({str(i): 1})
                else:
                    if str(i) not in result_dict[tttt_brcd]:
                        result_dict[tttt_brcd].update({str(i): 0})

    def get_dict_by_brcd_w_brcd_range_nSeq_from_FASTQ(self, fastq_list):
        result_dict = {}
        print("st : get_dict_by_brcd_w_brcd_range_nSeq_from_FASTQ")
        # self.brcd_obj is dict
        brcd_dict = self.brcd_obj

        for fastq_str in fastq_list:
            fastq_seq = fastq_str.upper()
            if self.rev_com_flag:
                fastq_seq = self.make_complement_string(fastq_str.upper())[::-1]

            en_brcd_i = self.brcd_pos_arr[1]
            if en_brcd_i < 0:
                en_brcd_i += len(fastq_seq)
                if en_brcd_i + self.brcd_len > len(fastq_seq):
                    en_brcd_i = len(fastq_seq) - self.brcd_len
            elif en_brcd_i == 0 or en_brcd_i + self.brcd_len > len(fastq_seq):
                en_brcd_i = len(fastq_seq) - self.brcd_len

            for brcd_i in range(self.brcd_pos_arr[0], en_brcd_i):
                candi_brcd = fastq_seq[brcd_i: brcd_i + self.brcd_len]

                if candi_brcd in brcd_dict:
                    self.check_barcode_freq_w_nSeq(brcd_dict[candi_brcd][0], fastq_seq, result_dict)

        print("DONE : get_dict_by_brcd_w_brcd_range_nSeq_from_FASTQ")
        return result_dict

    def get_dict_multi_p_seq_from_FASTQ(self, fastq_list):
        result_dict = {}
        print("get_dict_multi_p_seq_from_FASTQ starts ")
        # self.brcd_obj is list
        for val_arr in self.brcd_obj:
            tttt_brcd = val_arr[1].upper()
            ori_seq = val_arr[2].upper()
            edit_seq = val_arr[3].upper()
            for fastq_str in fastq_list:
                fastq_seq = fastq_str.upper()
                # check barcode
                if tttt_brcd in fastq_seq:
                    if tttt_brcd not in result_dict:
                        result_dict.update(
                            {tttt_brcd: {"Original sequence": 0, "Edited sequence": 0, "TTTT_Barcode_cnt": 1}})
                    else:
                        result_dict[tttt_brcd]["TTTT_Barcode_cnt"] += 1

                    # check seq after tttt_brcd seq
                    fastq_seq_aft_barcd = fastq_seq[fastq_seq.index(tttt_brcd) + len(tttt_brcd):]
                    # check original seq
                    if ori_seq in fastq_seq_aft_barcd:
                        result_dict[tttt_brcd]["Original sequence"] += 1

                    # check edited seq
                    if edit_seq in fastq_seq_aft_barcd:
                        result_dict[tttt_brcd]["Edited sequence"] += 1

        return result_dict

    def get_dict_multi_p_seq_from_whole_FASTQ(self, fastq_list):
        result_dict = {}
        print("get_dict_multi_p_seq_from_whole_FASTQ starts ")
        # self.brcd_obj is list
        for val_arr in self.brcd_obj:
            tttt_brcd = val_arr[1].upper()
            ori_seq = val_arr[2].upper()
            edit_seq = val_arr[3].upper()
            for fastq_str in fastq_list:
                fastq_seq = fastq_str.upper()
                # check barcode
                if tttt_brcd in fastq_seq:
                    if tttt_brcd not in result_dict:
                        result_dict.update(
                            {tttt_brcd: {"Original sequence": 0, "Edited sequence": 0, "TTTT_Barcode_cnt": 1}})
                    else:
                        result_dict[tttt_brcd]["TTTT_Barcode_cnt"] += 1

                    # check original seq
                    if ori_seq in fastq_seq:
                        result_dict[tttt_brcd]["Original sequence"] += 1

                    # check edited seq
                    if edit_seq in fastq_seq:
                        result_dict[tttt_brcd]["Edited sequence"] += 1

        return result_dict

    def get_dict_multi_p_4seq_from_whole_FASTQ(self, fastq_list):
        result_dict = {}
        print("get_dict_multi_p_4seq_from_whole_FASTQ starts ")
        # self.brcd_obj is list
        for val_arr in self.brcd_obj:
            tttt_brcd = val_arr[1].upper()
            wout_edit_seq = val_arr[2].upper()
            edit_seq = val_arr[3].upper()
            pos_1_seq = val_arr[4].upper()
            pos_2_seq = val_arr[5].upper()
            for fastq_str in fastq_list:
                fastq_seq = fastq_str.upper()
                # check barcode
                if tttt_brcd in fastq_seq:
                    if tttt_brcd not in result_dict:
                        result_dict.update(
                            {tttt_brcd: {"TTTT_Barcode_cnt": 1, "Target_sequences_without_edit": 0,
                                         "Target_sequences_with_edit_complete": 0,
                                         "Position_1_only": 0,
                                         "Position_2_only": 0}})
                    else:
                        result_dict[tttt_brcd]["TTTT_Barcode_cnt"] += 1

                    # check Target sequences without edit
                    if wout_edit_seq in fastq_seq:
                        result_dict[tttt_brcd]["Target_sequences_without_edit"] += 1

                    # check Target sequences with edit (complete)
                    if edit_seq in fastq_seq:
                        result_dict[tttt_brcd]["Target_sequences_with_edit_complete"] += 1

                    # check Position 1 only
                    if pos_1_seq in fastq_seq:
                        result_dict[tttt_brcd]["Position_1_only"] += 1

                    # check Position 2 only
                    if pos_2_seq in fastq_seq:
                        result_dict[tttt_brcd]["Position_2_only"] += 1

        return result_dict

    def merge_pool_list(self, pool_list):
        mege_dict = {}
        for data_dict in pool_list:
            for barcd_key, val_dict in data_dict.items():
                if barcd_key in mege_dict:
                    mege_dict[barcd_key]["TTTT_Barcode_cnt"] += val_dict["TTTT_Barcode_cnt"]
                    mege_dict[barcd_key]["Original sequence"] += val_dict["Original sequence"]
                    mege_dict[barcd_key]["Edited sequence"] += val_dict["Edited sequence"]
                else:
                    mege_dict.update({barcd_key: val_dict})

        return mege_dict

    def merge_4seq_pool_list(self, pool_list):
        mege_dict = {}
        for data_dict in pool_list:
            for barcd_key, val_dict in data_dict.items():
                if barcd_key in mege_dict:
                    mege_dict[barcd_key]["TTTT_Barcode_cnt"] += val_dict["TTTT_Barcode_cnt"]
                    mege_dict[barcd_key]["Target_sequences_without_edit"] += val_dict["Target_sequences_without_edit"]
                    mege_dict[barcd_key]["Target_sequences_with_edit_complete"] += val_dict["Target_sequences_with_edit_complete"]
                    mege_dict[barcd_key]["Position_1_only"] += val_dict["Position_1_only"]
                    mege_dict[barcd_key]["Position_2_only"] += val_dict["Position_2_only"]
                else:
                    mege_dict.update({barcd_key: val_dict})

        return mege_dict

    def get_dict_multi(self, bc_list):
        print("start get_dict_multi")
        # self.brcd_obj is list
        fastq_dir = self.brcd_obj[0]
        result_dict = {}
        for bc_arr in bc_list:
            fastq_nm = bc_arr[0]
            brcd_seq = bc_arr[1].upper()
            wt_seq = bc_arr[2].upper()
            edited_seq = bc_arr[3].upper()
            try:
                tmp = list(SeqIO.parse(fastq_dir + fastq_nm.replace(".gz", ""), "fastq"))
                seq_list = [[str(tmp[i].seq).upper(), str(tmp[i].reverse_complement().seq).upper()] for i in range(len(tmp))]

                for fastq_pair in seq_list:
                    fastq_seq = fastq_pair[0]
                    rvrs_com_fastq = fastq_pair[1]
                    if brcd_seq in fastq_seq:
                        if fastq_nm not in result_dict:
                            result_dict.update(
                                {fastq_nm: {"TTTT_Barcode_cnt": 1, "Target_sequences_without_edit": 0,
                                             "Target_sequences_with_edit_complete": 0,
                                             "Position_1_only": 0,
                                             "Position_2_only": 0}})
                        else:
                            result_dict[fastq_nm]["TTTT_Barcode_cnt"] += 1


                        if wt_seq in fastq_seq:
                            result_dict[fastq_nm]["Target_sequences_without_edit"] += 1
                        elif wt_seq in rvrs_com_fastq:
                            result_dict[fastq_nm]["Position_1_only"] += 1

                        if edited_seq in fastq_seq:
                            result_dict[fastq_nm]["Target_sequences_with_edit_complete"] += 1
                        elif edited_seq in rvrs_com_fastq:
                            result_dict[fastq_nm]["Position_2_only"] += 1
            except Exception as err:
                print("path : " + fastq_dir + fastq_nm)
                print(err)

        return result_dict

    def merge_pool_list_to_result_dict(self, pool_list, mege_dict):
        for data_dict in pool_list:
            for barcd_key, val_dict in data_dict.items():
                if barcd_key in mege_dict:
                    for key, val in val_dict.items():
                        mege_dict[barcd_key][key] += val
                else:
                    mege_dict.update({barcd_key: val_dict})

        # return mege_dict


