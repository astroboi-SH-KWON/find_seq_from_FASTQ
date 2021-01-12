import time
import multiprocessing as mp
import numpy as np
import os
import platform


import Util
import Logic
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
SYSTEM_NM = platform.system()
if SYSTEM_NM == 'Linux':
    # REAL
    pass
else:
    # DEV
    WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/"

FASTQ = 'FASTQ/'
BARCD_SEQ_FILE = "barcode_seq/210112_Novaseq_find_seq_from_FASTQ_index.txt"

MULTI_CNT = 10

############### end setting env ################

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
    print("process : " + str(mp.cpu_count()))
    pool = mp.Pool(processes=MULTI_CNT)
    ## analyze FASTQ seq after barcode seq
    # pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)
    ## analyze whole FASTQ seq
    pool_list = pool.map(logic.get_dict_multi_p_seq_from_whole_FASTQ, splited_fastq_list)

    merge_dict = logic.merge_pool_list(pool_list)

    util.make_dict_to_excel(WORK_DIR + "output/result_count", merge_dict)

def multi_processing_w_solo_fastq():
    util = Util.Utils()

    # fastq file name without ext
    solo_fastg_fl_nm_list = [
         "Monkey_PE_2K--TAATGCGC-GCCTCTAT_1"
        , "Monkey_PE_2K--GAATTCGT-AGGCTATA_1"
    ]

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)

    for fastq_fl_nm in solo_fastg_fl_nm_list:
        fastq_list = util.get_FASTQ_seq_list(WORK_DIR + FASTQ + fastq_fl_nm + '.fastq')

        # divide data_list by MULTI_CNT
        splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
        fastq_list.clear()

        print("process : " + str(mp.cpu_count()))
        pool = mp.Pool(processes=MULTI_CNT)
        ## analyze FASTQ seq after barcode seq
        # pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)
        ## analyze whole FASTQ seq
        pool_list = pool.map(logic.get_dict_multi_p_seq_from_whole_FASTQ, splited_fastq_list)

        merge_dict = logic.merge_pool_list(pool_list)
        pool.close()
        pool_list[:] = []

        util.make_dict_to_excel(WORK_DIR + "output/result_" + fastq_fl_nm, merge_dict)
        merge_dict.clear()

def multi_processing_test():
    util = Util.Utils()

    sources = util.get_files_from_dir(WORK_DIR + "FASTQ/" + '*.fastq')
    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)

    # key : D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq
    fastq_list = util.get_FASTQ_seq_dict(sources)['D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq']

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
    start_time = time.perf_counter()
    print("start >>>>>>>>>>>>>>>>>>")
    # multi_processing()
    multi_processing_w_solo_fastq()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))