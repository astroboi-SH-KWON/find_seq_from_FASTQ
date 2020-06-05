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