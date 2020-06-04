import pandas as pd
import Bio as bio
from Bio import SeqIO
from time import clock
import glob

import Util
import Logic
import LogicPrep
import Valid
############### start to set env ################
WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/"
# generates list of fastq files to analyze
SOURCES = glob.glob(WORK_DIR + "FASTQ/" + '*.fastq')

BARCD_SEQ_FILE = "barcode_seq/PE_library_input file_HW_GS_SJH_final_20200602.txt"


############### end setting env ################

def test():
    util = Util.Utils()

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)
    fastq_dict = util.get_FASTQ_seq(SOURCES)

    result_list = []
    for val_arr in brcd_list:
        for key, fasq_list in fastq_dict.items():
            for fastq_seq in fasq_list:
                if val_arr[1] in fastq_seq:
                    tmp_arr = [val_arr[0], val_arr[1]]
                    if val_arr[2] in fastq_seq:
                        tmp_arr.append(val_arr[2])
                    else:
                        tmp_arr.append("")
                    if val_arr[3] in fastq_seq:
                        tmp_arr.append(val_arr[3])
                    else:
                        tmp_arr.append("")
                    tmp_arr.append(fastq_seq)
                    result_list.append(tmp_arr)

    util.make_list_to_excel(WORK_DIR + "result", result_list)






start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
test()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))