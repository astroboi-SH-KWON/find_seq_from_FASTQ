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

def multi_processing_test_by_onefile():
    util = Util.Utils()

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)
    logic = Logic.Logics(brcd_list)

    fastq_list = util.get_FASTQ_seq(SOURCES)['D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq']

    # divide data_list by MULTI_CNT
    splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
    print("process : " + str(mp.cpu_count()))
    pool = mp.Pool(processes=MULTI_CNT)
    pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)

    util.make_dict_to_excel(WORK_DIR + "result_count", pool_list)

def multi_processing_test_by_onefile1():
    util = Util.Utils()

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)
    logic = Logic.Logics(brcd_list)

    fastq_list = util.get_FASTQ_seq(SOURCES)['D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq']

    # divide data_list by MULTI_CNT
    splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
    print("process : " + str(mp.cpu_count()))
    pool = mp.Pool(processes=MULTI_CNT)
    # result_list = pd.concat(pool.map(logic.get_multi_p_seq_from_FASTQ, splited_fastq_list))
    pool_list = pool.map(logic.get_list_multi_p_seq_from_FASTQ, splited_fastq_list)
    result_list = []
    for tmp_list in pool_list:
        result_list.extend(tmp_list)

    util.make_list_to_excel(WORK_DIR + "result_multi_processing", result_list)

if __name__ == '__main__':
    start_time = clock()
    print("start >>>>>>>>>>>>>>>>>>")
    multi_processing_test_by_onefile()
    print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))