import os
from pathlib import Path
import glob
import time
import json
import pickle

class dataSet():
    def __init__(self, disease):
        self._tpm = []
        self._fpkm = []
        self._fpkmuq = []
        self._name = disease
        self._average = []
        self._sd = []

    def addValue(self, fileds):
        self._tpm.append(float(fields[6]))
        self._fpkm.append(float(fields[7]))
        self._fpkmuq.append(float(fields[8]))

    @property
    def disease(self):
        return self.disease
    
    @property
    def tpm(self):
        return self._tpm
    
    @property
    def fpkm(self):
        return self._fpkm
    
    @property
    def fpkmuq(self):
        return self._fpkmuq
    



if __name__ == "__main__":
    start_time = time.time()

    # def project class
    GGT1 = dict()
    GGT2 = dict()
    GGT5 = dict()
    GGT6 = dict()
    GGT7 = dict()
    with open("tcga_abbr.txt", "r", encoding="utf-8") as k:
        while True:
            kline = k.readline()
            if not kline:
                break
            pj_name = kline.strip().split("\t")
            GGT1[pj_name[0]] = dataSet(pj_name[1])
            GGT2[pj_name[0]] = dataSet(pj_name[1])
            GGT5[pj_name[0]] = dataSet(pj_name[1])
            GGT6[pj_name[0]] = dataSet(pj_name[1])
            GGT7[pj_name[0]] = dataSet(pj_name[1])


    # get project id
    f = open("project_map.json", "r", encoding="utf-8")
    project = json.load(f)

    # get pwd
    rnaseq_file = "RNAseq"
    current_folder = os.path.join(os.getcwd(), rnaseq_file)
    #print("current file path:", current_folder)

    # get all file names
    files_and_folders = os.listdir(os.path.join(current_folder))
    print("all files:", files_and_folders)

    file_count = 0
    for filename in files_and_folders:
        tsv_files = glob.glob(os.path.join(current_folder, filename, "*.rna_seq.augmented_star_gene_counts.tsv"))
        if not tsv_files:
            continue
        count = 0
        #print("----------------------------------")
        with open(tsv_files[0], "r", encoding="utf-8") as file:
            file_count += 1
            while True:
                line = file.readline()  
                count += 1
                if not line:
                    #print("number of rows:", count)  
                    break
                if count == 2173:
                    fields = line.strip().split('\t')
                    GGT1[project[filename]].addValue(fields)
                if count == 6768:
                    fields = line.strip().split('\t')
                    GGT2[project[filename]].addValue(fields)
                if count == 2162:
                    fields = line.strip().split('\t')
                    GGT5[project[filename]].addValue(fields)
                if count == 12184:
                    fields = line.strip().split('\t')
                    GGT6[project[filename]].addValue(fields)
                if count == 6400:
                    fields = line.strip().split('\t')
                    GGT7[project[filename]].addValue(fields)
                    #print(fields[1]) 

    f.close()
    end_time = time.time()
    print("all time:", end_time - start_time)
    print("number of files:", file_count)

    # save results
    with open("result.pkl", "wb") as l:
        pickle.dump([GGT1, GGT2, GGT5, GGT6, GGT7], l)