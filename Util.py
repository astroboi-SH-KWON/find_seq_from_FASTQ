from os import listdir
from os.path import isfile, join
import pandas as pd
import openpyxl
from time import clock
import random
import math
import Bio as bio
from Bio import SeqIO
import glob

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    """
    reads the fastq files into a dictionary with the file names as keys
    :param
        FASTQ file
        @M00265:331:000000000-CJRDH:1:1101:19061:2012 1:N:0:200
        TGTTTCACAAACAAGACAGGGAACTGGCAGGCACCGAGGCATCTCTGCACCGAGGTGAAACAAGCTGCCATTTCATTACAGGCAAAGCTGAGCAAAAGTAGATATTACAAGACCAGCATGTACTCACCTCTCATGAAGCACTGTGGGTACGAAGGAAATGACTCAAATATGCTGTCTGAAGCCATCGCTTCCTCCTGAAAATGCACCCTCTTCTGAAGGCGGGGGACTCAATGATTTCTTTTACCTTCGGAGCGAAAACCAAGACAGGTCACTGTTTC
        +
        ABBBBFFFFFFAGCFFGFGGGCFFHHGHFGGFHHGGAECGDGHHHHGHFHHGCEGGHGHEHHHFHHHHHHHHHHBGHFHFFHGGEGHHHHEHHHHGHHHHHGHHHHFHHFHGHHHHHGHHHFHHHHHHHHHHHHHEHHHHHHHHHHGGHHGEAEGGHGHHGHHHHFFHHHFHHHHHHHEFGHHHHHGGHHGHHHHHHBFCBGHHHHHGHHGGGGGEB0FEGGAF;DDFFFFFFFFFFFFFFFFFFFFFFF?BFFFFFFFFFFFFFEF?FFFFFFFFFF
        @M00265:331:000000000-CJRDH:1:1101:18745:2023 1:N:0:200
        TCGGTCACAAACAAGACAGGGAACTGGCAGGCACCGAGGCATCTCTGCACCGAGGTGAAACAAGCTGCCATTTCATTACAGGCAAAGCTGAGCAAAAGTAGATATTACAAGACCAGCATGTACTCACCTCTCATGAAGCACTGTGGGTACGAAGGAAATGACTCAAATATGCTGTCTGAAGGCATCGCTTCCTCCTGAAAATGCACCCTCTTCTGAAGGCGGGGGACTCAATGATTTCTTTTACCTTCGGAGCGAAAACCAAGACAGGTCACTGTTTC
        +
        BABBBBBFFFFBGCBGFFGFGCFFHGBGFGGEHFGE?2EGCGGHHHFHFGHGGEGGGGHFHHEFGHHHH3FHHHFHHGG3BFFEAFGFGG3BFBBCGFHHHEGGGHGHFHBFHGHAG0GEGDGHHFHHHHHHBG121DGGHHCHFBAEFHDF/AA0FCGGDDGHCCGHHGFGHFHHHDDFDGCCGHADHHGHHHHHB:0:CCCGFGHGGGGGGGB0B09:A----BBFFFFFFFFFFFFFFFFFFFFFFDDBBBB=BBFFFE??FFFEFFFFFFFBFF
        @M00265:331:000000000-CJRDH:1:1101:19677:2030 1:N:0:200
        GACTTCACAAACAAGACAGGGAACTGGCAGGCACCGAGGCATCTCTGCACCGAGGTGAAACAAGCTGCCATTTCATTACAGGCAAAGCTGAGCAAAAGTAGATATTACAAGACCAGCATGTACTCACCTCTCATGAAGCACTGTGGGTACGAAGGAAATGACTCAAATATGCTGTCTGAAGCCATCGCTTCCTCCTGAAAATGCACCCTCTTCTGAAGGCGGGGGACTCAATGATTTCTTTTACCTTCGGAGCGAAAACCAAGACAGGTCACTGTTTC
        +

    :return
        {   'file name as key' : [FASTQ reads list from the file as values]
            , 'D:/000_WORK/KimHuiKwon/20200327/Fig 4d\\Fig4_d_RUNX1_6GC_rep3.fastq': 
    	        ['TGTTTCACAAACAAGACAGGGAACTGGCAGGCACCGAGGCATCTCTGCACCGAGGTGAAACAAGCTGCCATTTCATTACAGGCAAAGCTGAGCAAAAGTAGATATTACAAGACCAGCATGTACTCACCTCTCATGAAGCACTGTGGGTACGAAGGAAATGACTCAAATATGCTGTCTGAAGCCATCGCTTCCTCCTGAAAATGCACCCTCTTCTGAAGGCGGGGGACTCAATGATTTCTTTTACCTTCGGAGCGAAAACCAAGACAGGTCACTGTTTC'
    	        , 'TCGGTCACAAACAAGACAGGGAACTGGCAGGCACCGAGGCATCTCTGCACCGAGGTGAAACAAGCTGCCATTTCATTACAGGCAAAGCTGAGCAAAAGTAGATATTACAAGACCAGCATGTACTCACCTCTCATGAAGCACTGTGGGTACGAAGGAAATGACTCAAATATGCTGTCTGAAGGCATCGCTTCCTCCTGAAAATGCACCCTCTTCTGAAGGCGGGGGACTCAATGATTTCTTTTACCTTCGGAGCGAAAACCAAGACAGGTCACTGTTTC'
    	        , 'GACTTCACAAACAAGACAGGGAACTGGCAGGCACCGAGGCATCTCTGCACCGAGGTGAAACAAGCTGCCATTTCATTACAGGCAAAGCTGAGCAAAAGTAGATATTACAAGACCAGCATGTACTCACCTCTCATGAAGCACTGTGGGTACGAAGGAAATGACTCAAATATGCTGTCTGAAGCCATCGCTTCCTCCTGAAAATGCACCCTCTTCTGAAGGCGGGGGACTCAATGATTTCTTTTACCTTCGGAGCGAAAACCAAGACAGGTCACTGTTTC'
    	        , ... ]
    	    , 'D:/000_WORK/KimHuiKwon/20200327/Fig 4d\\Fig4_d_RUNX1_6GC_rep3.fastq':
    	        [...]
    	    ...}
    """
    def get_FASTQ_seq(self, sources):
        fastq_dict = {}
        for i in range(len(sources)):
            temp = list(SeqIO.parse(sources[i], "fastq"))
            fastq_dict[sources[i]] = [str(temp[k].seq) for k in range(len(temp))]

        return fastq_dict

    """
    :return
        result_list = [
                        [INDEX, TTTT_Barcode, Original sequence, Edited sequence]
                        , [CFTR_Forward TTT deletion, TTTTCTGCACGCATACTCTCAT, AAATATCATCTTTGGTGTTTCCTATGATGA, AAATATCATCGGTGTTTCCTATGATGAGTA]
                        ...
                    ]
    """
    def read_tb_txt(self, path):
        result_list = []
        with open(path, "r") as f:
            f.readline().replace("\n", "")
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == "":
                    break

                tmp_arr = tmp_line.split("\t")
                result_list.append(tmp_arr)

        return result_list

    """
        :return
            data_list = [
                            [INDEX, TTTT_Barcode, Original sequence, Edited sequence]
                            , [CFTR_Forward TTT deletion, TTTTCTGCACGCATACTCTCAT, AAATATCATCTTTGGTGTTTCCTATGATGA, ""]
                            , [CFTR_Forward TTT deletion, TTTTCTGCACGCATACTCTCAT, "", AAATATCATCGGTGTTTCCTATGATGAGTA]
                            , [CFTR_Forward TTT deletion, TTTTCTGCACGCATACTCTCAT, "", AAATATCATCGGTGTTTCCTATGATGAGTA]
                            ...
                        ]
        """
    def make_list_to_excel(self, result_path, data_list):
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        sheet.cell(row=row, column=1, value="index")
        sheet.cell(row=row, column=2, value='TTTT_Barcode')
        sheet.cell(row=row, column=3, value='Original sequence')
        sheet.cell(row=row, column=4, value='Edited sequence')
        sheet.cell(row=row, column=5, value='FASTQ sequence')

        for val_arr in data_list:
            row += 1
            col = 0
            for val in val_arr:
                col += 1
                sheet.cell(row=row, column=col, value=val)

        workbook.save(filename=result_path + self.ext_xlsx)

    def make_dict_to_excel(self, path, mege_dict):
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        sheet.cell(row=row, column=1, value="index")
        sheet.cell(row=row, column=2, value='TTTT_Barcode')
        sheet.cell(row=row, column=3, value='Original sequence')
        sheet.cell(row=row, column=4, value='Edited sequence')

        for barcd_key, val_dict in mege_dict.items():
            row += 1
            sheet.cell(row=row, column=1, value=str(row - 1))
            sheet.cell(row=row, column=2, value=barcd_key)
            sheet.cell(row=row, column=3, value=val_dict['Original sequence'])
            sheet.cell(row=row, column=4, value=val_dict['Edited sequence'])

        workbook.save(filename=path + self.ext_xlsx)

