from time import clock
import time
import glob
import multiprocessing as mp
import numpy as np
import os

import Util
import Logic
############### start to set env ################
# WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/"
WORK_DIR = os.getcwd() + "/"

BARCD_SEQ_FILE = "barcode_seq/PE_vivo_lib_input file_200609_final.txt"

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.9)

############### end setting env #################

def multi_processing():
    util = Util.Utils()

    sources = util.get_files_from_dir(WORK_DIR + "FASTQ/" + '*.fastq')
    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)

    fastq_list = []
    for sources_list in util.get_FASTQ_seq_dict(sources).values():
        fastq_list.extend(sources_list)

    # divide data_list by MULTI_CNT
    splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
    print("total cpu_count : " + str(TOTAL_CPU))
    print("will use : " + str(MULTI_CNT))
    pool = mp.Pool(processes=MULTI_CNT)
    # analyze FASTQ seq after barcode seq
    # pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)
    # analyze whole FASTQ seq
    pool_list = pool.map(logic.get_dict_multi_p_seq_from_whole_FASTQ, splited_fastq_list)

    merge_dict = logic.merge_pool_list(pool_list)

    util.make_dict_to_excel(WORK_DIR + "output/result_count", merge_dict)


if __name__ == '__main__':
    # start_time = clock()
    start_time = time.clock_gettime(1)
    print("start >>>>>>>>>>>>>>>>>>")
    multi_processing()
    # print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.clock_gettime(1) - start_time))
