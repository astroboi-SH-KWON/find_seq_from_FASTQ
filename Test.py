import pandas as pd
import Bio as bio
from Bio import SeqIO
from time import clock
import glob
import multiprocessing as mp
from threading import Thread
import numpy as np


import Util
import Logic
import LogicPrep
import Valid
############### start to set env ################
WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/"
# generates list of fastq files to analyze
SOURCES = glob.glob(WORK_DIR + "FASTQ/" + '*.fastq')

BARCD_SEQ_FILE = "barcode_seq/PE_library_input file_HW_GS_SJH_final_20200602.txt"

MULTI_CNT = 10


############### end setting env ################

def test():
    util = Util.Utils()
    logic = Logic.Logics()

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)
    fastq_dict = util.get_FASTQ_seq(SOURCES)

    result_list = logic.get_seq_from_FASTQ(brcd_list, fastq_dict)

    util.make_list_to_excel(WORK_DIR + "result", result_list)

# key : D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq
def multi_thread_test_by_onefile():
    util = Util.Utils()

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)
    logic = Logic.Logics(brcd_list)
    logic_prep = LogicPrep.LogicPreps()


    fastq_list = util.get_FASTQ_seq(SOURCES)['D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq']

    unit_len, remain = logic_prep.get_unit_len_n_remainer(fastq_list, MULTI_CNT)

    result_list = []
    # logic.get_multi_seq_from_FASTQ(brcd_list, fastq_list, [0, unit_len], result_list)

    th1 = Thread(target=logic.get_multi_seq_from_FASTQ, args=(fastq_list, [0, unit_len], result_list))
    th2 = Thread(target=logic.get_multi_seq_from_FASTQ, args=(fastq_list, [unit_len, unit_len * 2], result_list))
    th3 = Thread(target=logic.get_multi_seq_from_FASTQ, args=(fastq_list, [unit_len * 2, unit_len * 3], result_list))
    th4 = Thread(target=logic.get_multi_seq_from_FASTQ, args=(fastq_list, [unit_len * 3, unit_len * 4], result_list))
    th5 = Thread(target=logic.get_multi_seq_from_FASTQ, args=(fastq_list, [unit_len * 4, unit_len * 5 + remain], result_list))

    th1.start()
    th2.start()
    th3.start()
    th4.start()
    th5.start()

    th1.join()
    th2.join()
    th3.join()
    th4.join()
    th5.join()

    util.make_list_to_excel(WORK_DIR + "result_multi_thread", result_list)

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# test()
multi_thread_test_by_onefile()
# multi_processing_test_by_onefile()
# test2()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))