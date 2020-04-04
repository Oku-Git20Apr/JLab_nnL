#!/usr/bin/python3
####################################
## Make Missing Mass ##############
###################################



import sys
import time, os.path
from subprocess import call
#import concurrent.futures
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 
import numpy as np
import shlex
import subprocess
import re
import sys
import os,os.path
import glob
import datetime
nworkers=8


#####################################
date=datetime.date.today()
matrix="../matrix/matrix_new.list"
#matrix="../matrix/matrix_test.list"
#matrix="../matrix/matrix_wPathL.list"
root_dir="../rootfiles/mmass/ana_Lambda/"+str(date)
pdf_dir ="../pdf/mmass/ana_Lambda/"+str(date)
name=[]
param=[]
mode=[]
num=1
name.append("Lambda_small_OleH1")
mode.append("H1")
param.append("Lambda_H1.param")

name.append("nnL_small_Ole1")
mode.append("T")
param.append("nnL_1.param")

name.append("Lambda_small_OleH2")
mode.append("H1")
param.append("Lambda_H2.param")

name.append("Lambda_small_OleT")
mode.append("H2")
param.append("nnL_2.param")

name.append("nnL_small_Ole2")
mode.append("T")
param.append("nnL_2.param")

name.append("nnL_small_Ole3")
mode.append("T")
param.append("nnL_2.param")

name.append("nnL_small_Ole4")
mode.append("T")
param.append("nnL_3.param")

name.append("H3L_small_Ole")
mode.append("He")
param.append("nnL_3.param")

#param.append("f1_tuned_Lambda_twc.param")
#param.append("f1_Lambda_phase2_tuned.param")

#####################################


def Add():
    root_H=f'{root_dir}/Lambda_small_OleH_all.root'
    root_T=f'{root_dir}/nnL_small_Ole_all.root'
    add_H=f'hadd {root_H} {root_dir}/{name[0]}.root {root_dir}/{name[2]}.root'
    add_T=f'hadd {root_T} {root_dir}/{name[1]}.root {root_dir}/{name[4]}.root {root_dir}/{name[5]}.root {root_dir}/{name[6]}.root '
    print(add_H)
    subprocess.run([add_H],shell=True)
    print(add_T)
    subprocess.run([add_T],shell=True)
    
def ana(i):            
    names=name[i]
    MODE=mode[i]
    pfile=param[i]
    cmd=f'./bin/ana_Lambda -f ../run_list/nnlambda/{names}.list -p param/{pfile} -r {root_dir}/{names}.root -w {pdf_dir}/{names}.pdf -m {matrix} -{MODE}'
    print(cmd)
    subprocess.run([cmd],shell=True)

def main():

    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        for i in range(8):
            mappings = {executor.submit(ana,i)}
        
   
if os.path.exists(root_dir):
    while  num<100:
        if os.path.exists(root_dir+"_"+str(num)):
            num+=1
        else :
            root_dir=root_dir+"_"+str(num)
            pdf_dir=pdf_dir+"_"+str(num)
            break
        
print(f"Maked new root directory {root_dir}")
os.mkdir(root_dir)
print(f"Maked new pdf directory {pdf_dir}")
os.mkdir(pdf_dir)

main()
Add()

print()
print(f"=================< COMMENT >=================")
print(f"New Root directory : {root_dir}")
print(f"New pdf directory : {pdf_dir}")
print(f"==============================================")
