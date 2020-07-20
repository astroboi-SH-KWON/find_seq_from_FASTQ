import time
import glob
import multiprocessing as mp
import numpy as np
from Bio import SeqIO


import Util
import Logic
import LogicPrep
import Valid
############### start to set env ################
WORK_DIR = "D:/000_WORK/KimHuiKwon/20200720/WORK_DIR/"
# WORK_DIR = os.getcwd() + "/"

FASTQ_DIR = "David_NGS_PE2/"
BARCD_SEQ_FILE = "Reference_NGS.txt"
# BARCD_SEQ_FILE = "Reference_NGS_test.txt"

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.6)

############### end setting env ################

def multi_processing():
    util = Util.Utils()
    logic = Logic.Logics([WORK_DIR + FASTQ_DIR])

    brcd_list = util.read_tb_txt(WORK_DIR + BARCD_SEQ_FILE)
    # divide data_list by MULTI_CNT
    splited_bc_list = np.array_split(brcd_list, MULTI_CNT)
    print("total cpu_count : " + str(TOTAL_CPU))
    print("will use : " + str(MULTI_CNT))
    pool = mp.Pool(processes=MULTI_CNT)

    pool_list = pool.map(logic.get_dict_multi, splited_bc_list)
    merge_dict = logic.merge_4seq_pool_list(pool_list)
    util.make_4seq_dict_to_excel_(WORK_DIR + "output/result_count", merge_dict)

def test():
    # Fig2a_PE2_EMX1_9_rep1.fastq
    tmp = list(SeqIO.parse(WORK_DIR + FASTQ_DIR + "Fig2a_PE2_EMX1_9_rep1.fastq.gz".replace(".gz", ""), "fastq"))
    seq_list = [str(tmp[i].seq) for i in range(len(tmp))]
    print(seq_list)
    print(type(seq_list))



if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start >>>>>>>>>>>>>>>>>>")
    multi_processing()
    # test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))

