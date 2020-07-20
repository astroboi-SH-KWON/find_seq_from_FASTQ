import pandas as pd
import Bio as bio
from Bio import SeqIO
import time
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

BARCD_SEQ_FILE = "barcode_seq/PE_library_input file_HW_GS_SJH_final_20200602.txt"

MULTI_CNT = 10


############### end setting env ################

def test():
    util = Util.Utils()

    sources = util.get_files_from_dir(WORK_DIR + "FASTQ/" + '*.fastq')
    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)


    fastq_dict = util.get_FASTQ_seq_dict(sources)

    result_list = logic.get_seq_from_FASTQ(brcd_list, fastq_dict)

    util.make_list_to_excel(WORK_DIR + "result", result_list)

def multi_thread_test_by_onefile():
    util = Util.Utils()

    sources = util.get_files_from_dir(WORK_DIR + "FASTQ/" + '*.fastq')
    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)
    logic_prep = LogicPrep.LogicPreps()

    # key : D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq
    fastq_list = util.get_FASTQ_seq_dict(sources)['D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq']

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

def test2():
    result_count_dict = pd.read_excel(WORK_DIR + "result_count.xlsx").to_dict()['TTTT_Barcode']
    result_count_sv_dict = pd.read_excel(WORK_DIR + "result_count_from2100server_10p.xlsx").to_dict()['TTTT_Barcode']
    print("result_count_dict len : " + str(len(result_count_dict)))
    print("result_count_sv_dict len : " + str(len(result_count_sv_dict)))
    tmp_set1 = set()
    for val in result_count_dict.values():
        tmp_set1.add(val)
    print("result_count_dict len : " + str(len(tmp_set1)))
    tmp_set2 = set()
    for val2 in result_count_sv_dict.values():
        tmp_set2.add(val2)
    print("result_count_sv_dict len : " + str(len(tmp_set2)))

def test3():
    num_arr = [0.916, 0.978, 0.658, 0.794, 0.995, 0.955, 0.818]
    print(sorted(num_arr))
    print(str((len(num_arr) + 1) / 4))


start_time = time.perf_counter()
print("start >>>>>>>>>>>>>>>>>>")
# test()
# multi_thread_test_by_onefile()
test3()
print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))