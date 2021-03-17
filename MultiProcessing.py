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
    # WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/"
    WORK_DIR = "D:/000_WORK/LeeMinYung/20210316/WORK_DIR/"  # 20210316

FASTQ = 'FASTQ/'
# BARCD_SEQ_FILE = "barcode_seq/210112_Novaseq_find_seq_from_FASTQ_index.txt"
# BARCD_SEQ_FILE = "barcode_seq/210113_Novaseq_find_seq_from_FASTQ_index.txt"
# BARCD_SEQ_FILE = "barcode_seq/210316_find seq_MY_pig_PE_1st analysis.txt"  # 20210316
BARCD_SEQ_FILE = "barcode_seq/210316_find seq_MY_pig_PE_2nd analysis_pam coedited.txt"  # 20210316
# for multi_processing_w_brcd_arr()
BARCD_FL_ARR = ["barcode_seq/210317_find seq_MY_1st.txt", "barcode_seq/210317_find seq_MY_2nd.txt"]  # 20210317

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.5)

############### end setting env ################

def multi_processing():
    util = Util.Utils()

    print("read FASTQ files path")
    sources = util.get_files_from_dir(WORK_DIR + "FASTQ/" + '*.fastq')
    print("read barcode list")
    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)

    print("get_FASTQ_seq_dict")
    fastq_list = []
    for sources_list in util.get_FASTQ_seq_dict(sources).values():
        fastq_list.extend(sources_list)

    # divide data_list by MULTI_CNT
    splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
    print("platform.system() : ", SYSTEM_NM)
    print("total cpu_count : ", str(TOTAL_CPU))
    print("will use : ", str(MULTI_CNT))
    pool = mp.Pool(processes=MULTI_CNT)
    ## analyze FASTQ seq after barcode seq
    pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)
    ## analyze whole FASTQ seq
    # pool_list = pool.map(logic.get_dict_multi_p_seq_from_whole_FASTQ, splited_fastq_list)

    merge_dict = logic.merge_pool_list(pool_list)

    util.make_dict_to_excel(
        WORK_DIR + "output/result_count_" + BARCD_SEQ_FILE.replace("barcode_seq/", "").replace(".txt", ""), merge_dict)


def multi_processing_w_brcd_arr():  # 20210317
    util = Util.Utils()
    fastq_fl_arr = [
                    "503_707.fastq"
                    , "504_708.fastq"
                    ]

    for brcd_fl in BARCD_FL_ARR:
        print("read barcode list")
        brcd_list = util.read_tb_txt(WORK_DIR + brcd_fl)

        logic = Logic.Logics(brcd_list)

        for fastq_fl in fastq_fl_arr:
            print("get_FASTQ_seq_list : ", fastq_fl)
            fastq_list = util.get_FASTQ_seq_list(WORK_DIR + FASTQ + fastq_fl)

            # divide data_list by MULTI_CNT
            splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
            print("platform.system() : ", SYSTEM_NM)
            print("total cpu_count : ", str(TOTAL_CPU))
            print("will use : ", str(MULTI_CNT))
            pool = mp.Pool(processes=MULTI_CNT)
            ## analyze FASTQ seq after barcode seq
            pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)
            ## analyze whole FASTQ seq
            # pool_list = pool.map(logic.get_dict_multi_p_seq_from_whole_FASTQ, splited_fastq_list)

            merge_dict = logic.merge_pool_list(pool_list)
            pool.close()
            pool_list[:] = []

            util.make_dict_to_excel(
                WORK_DIR + "output/result_count_" + fastq_fl.replace(".fastq", "") + "_" + brcd_fl.replace(
                    "barcode_seq/", "").replace(".txt", ""), merge_dict)


def multi_processing_w_solo_fastq():
    util = Util.Utils()

    # fastq file name without ext
    solo_fastq_fl_nm_list = [
         "Monkey_PE_2K_rep2_add"
        , "Monkey_PE_2K_Rep1_add"
        , "Monkey_PE_2K_Un"
    ]

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)

    for fastq_fl_nm in solo_fastq_fl_nm_list:
        fastq_list = util.get_FASTQ_seq_list(WORK_DIR + FASTQ + fastq_fl_nm + '.fastq')

        # divide data_list by MULTI_CNT
        splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
        fastq_list.clear()

        print("platform.system() : ", SYSTEM_NM)
        print("total cpu_count : ", str(TOTAL_CPU))
        print("will use : ", str(MULTI_CNT))
        pool = mp.Pool(processes=MULTI_CNT)
        ## analyze FASTQ seq after barcode seq
        pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)
        ## analyze whole FASTQ seq
        # pool_list = pool.map(logic.get_dict_multi_p_seq_from_whole_FASTQ, splited_fastq_list)

        merge_dict = logic.merge_pool_list(pool_list)
        pool.close()
        pool_list[:] = []

        util.make_dict_to_excel(WORK_DIR + "output/result_" + fastq_fl_nm, merge_dict)
        merge_dict.clear()


def multi_processing_split_big_files_then_find_seq_from_FASTQ():
    print('multi_processing_split_big_files_then_find_seq_from_FASTQ')
    util = Util.Utils()

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)
    logic = Logic.Logics(brcd_list)

    # fastq file name without ext
    big_fastq_fl_nm_list = [
        "MY_pig_PE_NG"
        , "MY_pig_PE_NGG"
    ]
    fastq_ext = '.fastq'

    for fastq_fl_nm in big_fastq_fl_nm_list:
        # split big file
        split_init = {'big_file_path': WORK_DIR + FASTQ + fastq_fl_nm + fastq_ext
                            , 'num_row': 4000000
                            , 'splited_files_dir': WORK_DIR + FASTQ + fastq_fl_nm + "/"
                            , 'output_file_nm': fastq_fl_nm
                            , 'output_file_ext': fastq_ext

        }
        util.split_big_file_by_row(split_init)

        # get splited_files path
        sources = util.get_files_from_dir(split_init['splited_files_dir'] + '*.fastq')

        result_dict = {}
        for splited_fastq_fl in sources:
            print("get_FASTQ_seq_list :", splited_fastq_fl)
            fastq_list = util.get_FASTQ_seq_list(splited_fastq_fl)

            # divide data_list by MULTI_CNT
            splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
            fastq_list.clear()

            print("platform.system() : ", SYSTEM_NM)
            print("total cpu_count : ", str(TOTAL_CPU))
            print("will use : ", str(MULTI_CNT))
            pool = mp.Pool(processes=MULTI_CNT)
            ## analyze FASTQ seq after barcode seq
            pool_list = pool.map(logic.get_dict_multi_p_seq_from_FASTQ, splited_fastq_list)
            ## analyze whole FASTQ seq
            # pool_list = pool.map(logic.get_dict_multi_p_seq_from_whole_FASTQ, splited_fastq_list)

            print("merge_pool_list_to_result_dict")
            logic.merge_pool_list_to_result_dict(pool_list, result_dict)
            pool.close()
            pool_list[:] = []

        print("make excel file")
        util.make_dict_to_excel(
            WORK_DIR + "output/result_" + fastq_fl_nm + "_" + BARCD_SEQ_FILE.replace("barcode_seq/", "").replace(".txt",
                                                                                                                 ""),
            result_dict)
        result_dict.clear()


def multi_processing_test():
    util = Util.Utils()

    sources = util.get_files_from_dir(WORK_DIR + "FASTQ/" + '*.fastq')
    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)

    logic = Logic.Logics(brcd_list)

    # key : D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq
    fastq_list = util.get_FASTQ_seq_dict(sources)['D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/FASTQ\\18.fastq']

    # divide data_list by MULTI_CNT
    splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
    print("platform.system() : ", SYSTEM_NM)
    print("total cpu_count : ", str(TOTAL_CPU))
    print("will use : ", str(MULTI_CNT))
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
    # multi_processing_w_solo_fastq()
    # multi_processing_split_big_files_then_find_seq_from_FASTQ()
    multi_processing_w_brcd_arr()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))