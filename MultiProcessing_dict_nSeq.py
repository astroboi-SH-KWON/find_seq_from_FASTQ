import time
import multiprocessing as mp
import numpy as np
import os
import platform


import Util
import Logic
import LogicPrep
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

FASTQ = '/FASTQ/'
BRCD = '/barcode_seq/'
# for multi_processing_w_brcd_arr()
BARCD_FL_ARR = [
                # "EB_BE_findseq.txt"
                # "PE_vivo_lib_input file_200609_final.txt"
                "1k_Barcode_silent.txt"
                # , "2104222.txt"
                # , "2104223.txt"
                # , "2104225.txt"
                # , "2104227.txt"
                # , "2104221.txt"
                ]

# fastq file name without ext
FQ_EXT = '.fastq'
# FQ_EXT = '.fq'
FASTG_ARR = [
            # "ABEng_1"
            "ABEngg"
             , "ABEng"
            ]

BRCD_POS_ARR_fr_FQ = [0, 0]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.9)

############### end setting env ################


# multi BIG fastq files with multi barcode files
def multi_processing_dict_split_big_files_then_find_seq_from_FASTQ_w_brcd_arr():
    print('st multi_processing_dict_split_big_files_then_find_seq_from_FASTQ_w_brcd_arr')
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    for brcd_path in BARCD_FL_ARR:
        brcd_list = util.read_tb_txt(WORK_DIR + BRCD + brcd_path)
        brcd_dict = logic_prep.list_to_dict_by_ele_key(brcd_list, 1)
        logic = Logic.Logics(brcd_dict)

        ### st set logic init ####################################
        # set bar_code estimated start and end position from fastq, default : [0, 0]
        logic.brcd_pos_arr = BRCD_POS_ARR_fr_FQ
        # set bar_code length, default : 20
        logic.brcd_len = logic_prep.check_brcd_length(brcd_dict, brcd_path)
        print("BRCD_len :", logic.brcd_len)
        # set reversed complement flag
        # logic.rev_com_flag = True
        ### en set logic init ####################################

        # BIG fastq files
        for fastq_fl_nm in FASTG_ARR:
            # # st split_big_file_by_row
            split_init = {'big_file_path': WORK_DIR + FASTQ + fastq_fl_nm + FQ_EXT
                                , 'num_row': 4000000
                                , 'splited_files_dir': WORK_DIR + FASTQ + fastq_fl_nm + "/"
                                , 'output_file_nm': fastq_fl_nm
                                , 'output_file_ext': FQ_EXT

            }
            # util.split_big_file_by_row(split_init)
            # # DONE split_big_file_by_row

            # get splited_files path
            sources = util.get_files_from_dir(split_init['splited_files_dir'] + '*' + FQ_EXT)
            # check sources' length
            if len(sources) == 0:
                print("[ERROR] : check *." + FQ_EXT, "files path, it's empty")
                exit()

            result_dict = {}
            for splited_fastq_fl in sources:
                fastq_list = util.get_FASTQ_seq_list(splited_fastq_fl)

                # divide data_list by MULTI_CNT
                splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
                fastq_list.clear()

                print("platform.system() : ", SYSTEM_NM)
                print("total cpu_count : ", str(TOTAL_CPU))
                print("will use : ", str(MULTI_CNT))
                pool = mp.Pool(processes=MULTI_CNT)
                ## analyze FASTQ seq after barcode seq
                pool_list = pool.map(logic.get_dict_by_brcd_w_brcd_range_from_FASTQ, splited_fastq_list)

                print("merge_pool_list_to_result_dict")
                logic.merge_pool_list_to_result_dict(pool_list, result_dict)
                pool.close()
                pool_list[:] = []

            logic_prep.add_missing_brcd_to_dict(brcd_list, result_dict)
            print("make excel file")
            util.make_dict_to_excel(
                WORK_DIR + "output/result_" + fastq_fl_nm + "_" + brcd_path.replace("barcode_seq/", "").replace(".txt",
                                                                                                                     ""),
                result_dict)
            result_dict.clear()


def multi_processing_dict_w_brcd_arr_nSeq():
    print('st multi_processing_dict_w_brcd_arr_nSeq')
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    for brcd_fl in BARCD_FL_ARR:
        print("read barcode list")
        brcd_list = util.read_tb_txt(WORK_DIR + BRCD + brcd_fl)
        brcd_dict = logic_prep.list_to_dict_by_ele_key(brcd_list, 1)
        logic = Logic.Logics(brcd_dict)

        ### st set logic init ####################################
        # # set bar_code estimated start and end position from fastq, default : [0, 0]
        logic.brcd_pos_arr = BRCD_POS_ARR_fr_FQ
        # # set bar_code length, default : 20
        logic.brcd_len = logic_prep.check_brcd_length(brcd_dict, brcd_fl)
        print("BRCD_len :", logic.brcd_len)
        # # set reversed complement flag
        # logic.rev_com_flag = True
        ### en set logic init ####################################

        for fastq_fl in FASTG_ARR:
            result_dict = {}
            print("get_FASTQ_seq_list : ", fastq_fl)
            fastq_list = util.get_FASTQ_seq_list(WORK_DIR + FASTQ + fastq_fl + FQ_EXT)

            # divide data_list by MULTI_CNT
            splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
            print("platform.system() : ", SYSTEM_NM)
            print("total cpu_count : ", str(TOTAL_CPU))
            print("will use : ", str(MULTI_CNT))
            pool = mp.Pool(processes=MULTI_CNT)
            ## analyze FASTQ seq after barcode seq
            pool_list = pool.map(logic.get_dict_by_brcd_w_brcd_range_nSeq_from_FASTQ, splited_fastq_list)

            print("merge_pool_list_to_result_dict")
            logic.merge_pool_list_to_result_dict(pool_list, result_dict)
            pool.close()
            pool_list[:] = []

            # logic_prep.add_missing_brcd_to_dict(brcd_list, result_dict)
            util.make_dict_to_excel_nSeq(WORK_DIR + "output/result_count_" + fastq_fl + "_" + brcd_fl.replace(".txt", ""),
                                    result_dict)

def split_big_file():
    print('st split_big_file')
    util = Util.Utils()
    for fastq_fl_nm in FASTG_ARR:
        # split big file
        split_init = {'big_file_path': WORK_DIR + FASTQ + fastq_fl_nm + FQ_EXT
            , 'num_row': 4000000
            , 'splited_files_dir': WORK_DIR + FASTQ + fastq_fl_nm + "/"
            , 'output_file_nm': fastq_fl_nm
            , 'output_file_ext': FQ_EXT

                      }
        util.split_big_file_by_row(split_init)
    print('DONE split_big_file')


def test():
    brcd_seq = "ABCDEFG"

    for i in range(2, len(brcd_seq)):
        print(brcd_seq[i])


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start >>>>>>>>>>>>>>>>>>")
    # test()
    multi_processing_dict_w_brcd_arr_nSeq()
    # multi_processing_dict_split_big_files_then_find_seq_from_FASTQ_w_brcd_arr()
    # split_big_file()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))