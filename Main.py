import time


import Util
import Logic
import LogicPrep
import Valid
############### start to set env ################
# WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/"
WORK_DIR = "D:/000_WORK/KimHuiKwon/20200618/WORK_DIR/"

SUB_OUT_DIR = "total/"
############### end setting env #################



def merge_multi_processing_excel_result():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    txt_sources = util.get_files_from_dir(WORK_DIR + "output/" + '*.txt')

    total_list = []
    for txt_file in txt_sources:
        total_list.extend(util.read_tb_txt(txt_file))

    merge_dict = logic_prep.make_list_to_dict(total_list)

    util.make_dict_to_excel(WORK_DIR + "output/merge_result_count", merge_dict)

def merge_multi_processing_4seq_excel_result():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    txt_sources = util.get_files_from_dir(WORK_DIR + "output/" + SUB_OUT_DIR + '*.txt')

    total_list = []
    for txt_file in txt_sources:
        total_list.extend(util.read_tb_txt(txt_file))

    merge_dict = logic_prep.make_4seq_list_to_dict(total_list)

    util.make_4seq_dict_to_excel(WORK_DIR + "output/" + SUB_OUT_DIR + "merge_result_count", merge_dict)


start_time = time.perf_counter()
print("start >>>>>>>>>>>>>>>>>>")
# merge_multi_processing_excel_result()
merge_multi_processing_4seq_excel_result()
# merge_multi_processing_4seq_excel_result()
print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))