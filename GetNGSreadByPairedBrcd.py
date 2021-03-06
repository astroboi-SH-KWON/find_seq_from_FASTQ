import time
import os
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform
import operator

import Util
#################### st env ####################
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
SYSTEM_NM = platform.system()

if SYSTEM_NM == 'Linux':
    # REAL
    FASTQ_DIR = '/media/backup/minyoung_precancer_1st_backup/run_CRISPResso2_in_cmd_1st/'
else:
    # DEV
    WORK_DIR = "D:/000_WORK/JangHyeWon_LeeMinYung/20200703/WORK_DIR/"
    FASTQ_DIR = "D:/000_WORK/JangHyeWon_LeeMinYung/20200703/WORK_DIR/"

FASTQ = 'FASTQ/'
IN = 'input/'
OU = 'output/'

MIS_MATCH_CUT_OFF = 1

TRGT_FASTQ_NM = '253_1_S6_L001_R1_001_join'
TRGT_FASTQ_EXT = '.fastq'
BRCD_FILE = "barcode_list.txt"
PAIRED_BRCD_FILE = "paied_brcd_20.txt"
BRCD1_POS =[0, 9]
LEN_BRCD = 9
GAP_ARR = [0, 1, 2, 3]
CONST1 = "AGTACGTACGAGTC"  # 14 bp
CONST2 = "GTACTCGCAGTAGTC"  # 15 bp
CONST2_POS = [30, 37]

os.makedirs(WORK_DIR + IN, exist_ok=True)
os.makedirs(WORK_DIR + OU, exist_ok=True)

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)

brcd_list = Util.Utils().read_tsv_ignore_N_line(WORK_DIR + IN + BRCD_FILE)
BRCD_ARR = [tmp_arr[0] for tmp_arr in brcd_list]

paied_brcd_list = Util.Utils().read_tsv_ignore_N_line(WORK_DIR + IN + PAIRED_BRCD_FILE)
PAIRED_BRCD_ARR = [tmp_arr[0] + '_' + tmp_arr[1] for tmp_arr in paied_brcd_list]
#################### en env ####################


def multi_NGS_read_by_PairedBrcd_aftr_split():
    print("multi_NGS_read_by_PairedBrcd_aftr_split", PAIRED_BRCD_FILE)
    util = Util.Utils()

    fastq_dir = TRGT_FASTQ_NM + "/"
    sources = util.get_files_from_dir(FASTQ_DIR + FASTQ + fastq_dir + TRGT_FASTQ_NM + '*' + TRGT_FASTQ_EXT)

    result_list = []
    for fastq_path in sources:
        temp = list(SeqIO.parse(fastq_path, "fastq"))
        fastq_list = [str(temp[k].seq) for k in range(len(temp))]
        temp.clear()

        # divide data_list by MULTI_CNT
        splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
        fastq_list.clear()
        print("platform.system() : ", SYSTEM_NM)
        print("total cpu_count : ", str(TOTAL_CPU))
        print("will use : ", str(MULTI_CNT))
        pool = mp.Pool(processes=MULTI_CNT)

        pool_list = pool.map(get_NGS_read_by_selected_PairedBrcd, splited_fastq_list)
        pool.close()

        # merge results
        print('\n\nmerge results', fastq_path)
        for tmp_list in pool_list:
            result_list.extend(tmp_list)
        pool_list.clear()

    # sort by paired_brcd
    result_list.sort(key=operator.itemgetter(0))

    result_fl_path = FASTQ_DIR + FASTQ + fastq_dir + 'tot_result_' + TRGT_FASTQ_NM + '.txt'
    with open(result_fl_path, 'w') as f_ou:
        print('make result file', result_fl_path)
        for result_arr in result_list:
            tmp_line = ''
            for tmp_val in result_arr:
                tmp_line = tmp_line + tmp_val + '\t'
            f_ou.write(tmp_line[:-1] + '\n')
    result_list.clear()


def check_brcd_by_brcd1(brcd1_fr_ngs):
    for paired_brcd in PAIRED_BRCD_ARR:
        if brcd1_fr_ngs == paired_brcd[0: len(brcd1_fr_ngs)]:
            return paired_brcd
    return ''


def get_NGS_read_by_selected_PairedBrcd(fastq_list):
    print("get_NGS_read_by_selected_PairedBrcd", PAIRED_BRCD_FILE)
    result_list = []
    for ngs_read in fastq_list:
        brcd1_fr_ngs = ngs_read[BRCD1_POS[0]: BRCD1_POS[1]]
        paired_brcd = check_brcd_by_brcd1(brcd1_fr_ngs)
        if paired_brcd == '':
            continue

        cons2_pos = ngs_read.find(CONST2, BRCD1_POS[1])
        len_cons2 = len(CONST2)

        if CONST2_POS[0] <= cons2_pos <= CONST2_POS[1]:
            brcd2_fr_ngs = ngs_read[cons2_pos - LEN_BRCD: cons2_pos]
            brcd_paired_key = brcd1_fr_ngs + "_" + brcd2_fr_ngs
            if paired_brcd != brcd_paired_key:
                continue

            result_list.append([brcd_paired_key, ngs_read])

        else:  # get n_mismatch const2 by MIS_MATCH_CUT_OFF
            for i in range(CONST2_POS[0], CONST2_POS[1]):
                const_candi = ngs_read[i: i + len_cons2].upper()

                # check length
                if len(const_candi) != len_cons2:
                    continue

                cnt = 0
                for idx in range(len_cons2):
                    if const_candi[idx] != CONST2[idx]:
                        cnt += 1
                    if cnt > MIS_MATCH_CUT_OFF:
                        break
                if cnt <= MIS_MATCH_CUT_OFF:
                    brcd2_fr_ngs = ngs_read[i - LEN_BRCD: i]
                    brcd_paired_key = brcd1_fr_ngs + "_" + brcd2_fr_ngs
                    if paired_brcd != brcd_paired_key:
                        continue

                    result_list.append([brcd_paired_key, ngs_read])
    return result_list


def get_NGS_read_by_1500x1500_cell_id(fastq_list):
    print('get_NGS_read_by_1500x1500_cell_id')
    result_list = []
    for ngs_read in fastq_list:
        brcd1_fr_ngs = ngs_read[BRCD1_POS[0]: BRCD1_POS[1]]
        if brcd1_fr_ngs not in BRCD_ARR:
            continue

        cons2_pos = ngs_read.find(CONST2, BRCD1_POS[1])
        len_cons2 = len(CONST2)

        if CONST2_POS[0] <= cons2_pos <= CONST2_POS[1]:
            brcd2_fr_ngs = ngs_read[cons2_pos - LEN_BRCD: cons2_pos]
            if brcd2_fr_ngs not in BRCD_ARR:
                continue

            brcd_paired_key = brcd1_fr_ngs + "_" + brcd2_fr_ngs
            result_list.append([brcd_paired_key, ngs_read])

        else:  # get n_mismatch const2 by MIS_MATCH_CUT_OFF
            for i in range(CONST2_POS[0], CONST2_POS[1]):
                const_candi = ngs_read[i: i + len_cons2].upper()

                # check length
                if len(const_candi) != len_cons2:
                    continue

                cnt = 0
                for idx in range(len_cons2):
                    if const_candi[idx] != CONST2[idx]:
                        cnt += 1
                    if cnt > MIS_MATCH_CUT_OFF:
                        break
                if cnt <= MIS_MATCH_CUT_OFF:
                    brcd2_fr_ngs = ngs_read[i - LEN_BRCD: i]
                    if brcd2_fr_ngs not in BRCD_ARR:
                        continue

                    brcd_paired_key = brcd1_fr_ngs + "_" + brcd2_fr_ngs
                    result_list.append([brcd_paired_key, ngs_read])
    return result_list


def multi_NGS_read_by_brcd_aftr_split():
    print('multi_NGS_read_by_brcd_aftr_split\n')
    util = Util.Utils()

    fastq_dir = TRGT_FASTQ_NM + "/"
    sources = util.get_files_from_dir(FASTQ_DIR + FASTQ + fastq_dir + TRGT_FASTQ_NM + '*' + TRGT_FASTQ_EXT)

    for fastq_path in sources:
        temp = list(SeqIO.parse(fastq_path, "fastq"))
        fastq_list = [str(temp[k].seq) for k in range(len(temp))]
        temp.clear()

        # divide data_list by MULTI_CNT
        splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
        fastq_list.clear()
        print("platform.system() : ", SYSTEM_NM)
        print("total cpu_count : ", str(TOTAL_CPU))
        print("will use : ", str(MULTI_CNT))
        pool = mp.Pool(processes=MULTI_CNT)

        pool_list = pool.map(get_NGS_read_by_1500x1500_cell_id, splited_fastq_list)
        pool.close()

        # merge results
        print('\n\nmerge results', fastq_path)
        result_list = []
        for tmp_list in pool_list:
            result_list.extend(tmp_list)
        pool_list.clear()

        # sort by paired_brcd
        result_list.sort(key=operator.itemgetter(0))

        result_fl_path = fastq_path.replace(TRGT_FASTQ_EXT, '_result.txt')
        with open(result_fl_path, 'w') as f_ou:
            print('make result file', result_fl_path)
            for result_arr in result_list:
                tmp_line = ''
                for tmp_val in result_arr:
                    tmp_line = tmp_line + tmp_val + '\t'
                f_ou.write(tmp_line[:-1] + '\n')
        result_list.clear()


def split_trgt_fastq_file():
    print(split_trgt_fastq_file)
    util = Util.Utils()
    init_split_file = {'big_file_path': FASTQ_DIR + FASTQ + TRGT_FASTQ_NM + TRGT_FASTQ_EXT
        , 'num_row': 4000000
        , 'splited_files_dir': FASTQ_DIR + FASTQ + TRGT_FASTQ_NM + "/"
        , 'output_file_nm': TRGT_FASTQ_NM
        , 'output_file_ext': TRGT_FASTQ_EXT
                       }
    util.split_big_file_by_row(init_split_file)


def multi_make_files_by_brcd_from_NGS_read_by_brcd_aftr_split():
    paired_brcd_arr = [brcd1 + '_' + brcd2 for brcd2 in BRCD_ARR for brcd1 in BRCD_ARR]

    # divide data_list by MULTI_CNT
    splited_paired_brcd_arr = np.array_split(paired_brcd_arr, MULTI_CNT)
    paired_brcd_arr.clear()
    print("platform.system() : ", SYSTEM_NM)
    print("total cpu_count : ", str(TOTAL_CPU))
    print("will use : ", str(MULTI_CNT))
    pool = mp.Pool(processes=MULTI_CNT)

    pool.map(make_files_by_brcd, splited_paired_brcd_arr)
    pool.close()


def make_files_by_brcd(paired_brcd_arr):
    util = Util.Utils()
    fastq_dir = TRGT_FASTQ_NM + "/"
    sources = util.get_files_from_dir(FASTQ_DIR + FASTQ + fastq_dir + TRGT_FASTQ_NM + '*_result.txt')

    for paired_brcd in paired_brcd_arr:
        result_list = []
        for result_fl_path in sources:
            brcd_ngs_read_list = util.read_tsv_ignore_N_line(result_fl_path, 0)
            for brcd_ngs_read_arr in brcd_ngs_read_list:
                if paired_brcd != brcd_ngs_read_arr[0]:
                    continue
                result_list.append(brcd_ngs_read_arr)

        if len(result_list) == 0:
            continue

        result_dir_path = FASTQ_DIR + FASTQ + fastq_dir + 'result_by_1500x1500_brcd/'
        os.makedirs(result_dir_path, exist_ok=True)
        with open(result_dir_path + paired_brcd + '.txt', 'w') as f_ou:
            print('make file for paired barcode [', paired_brcd, ']')
            for result_arr in result_list:
                tmp_line = ''
                for tmp_val in result_arr:
                    tmp_line = tmp_line + tmp_val + '\t'
                f_ou.write(tmp_line[:-1] + '\n')
        result_list.clear()


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # split_trgt_fastq_file()
    # multi_NGS_read_by_brcd_aftr_split()
    # multi_make_files_by_brcd_from_NGS_read_by_brcd_aftr_split()
    multi_NGS_read_by_PairedBrcd_aftr_split()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))