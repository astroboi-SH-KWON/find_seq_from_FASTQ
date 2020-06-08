from time import clock


import Util
import Logic
import LogicPrep
import Valid
############### start to set env ################
WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200604/WORK_DIR/"



def merge_multi_processing_excel_result():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    excel_sources = util.get_files_from_dir(WORK_DIR + "output/" + '*.txt')

    total_list = []
    for txt_file in excel_sources:
        total_list.extend(util.read_tb_txt(txt_file))

    merge_dict = logic_prep.make_list_to_dict(total_list)

    util.make_dict_to_excel(WORK_DIR + "output/merge_result_count", merge_dict)


start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
merge_multi_processing_excel_result()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))