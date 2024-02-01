#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 15:17:35 2021

@author: serena-mo
"""

import numpy as np
import random
from collections import Counter

class Error:
    def __init__(self, strands, pr, sub_pr, del_pr, ins_pr):
        """
        Initialize the Error class.

        Parameters:
        - strands: List of input sequences (each sequence is a NumPy array).
        - pr: Probability of dropout (sequence not sampled).
        - sub_pr: Probability of substitution error.
        - ins_pr: Probability of insertion error.
        - del_pr: Probability of deletion error.
        """
        self.strands = strands
        self.pr = pr
        self.sub_pr = sub_pr
        self.ins_pr = ins_pr
        self.del_pr = del_pr

    def random_error(self, y):
        """
        Introduce random errors into a sequence.

        Parameters:
        - y: Input sequence (NumPy array).

        Returns:
        - af_code_array: Sequence after introducing errors (NumPy array).
        - del_num: Number of deletion errors.
        - ins_num: Number of insertion errors.
        - sub_num: Number of substitution errors.
        """
        y = np.copy(y)
        sub_pr = self.sub_pr
        ins_pr = self.ins_pr
        del_pr = self.del_pr
        pt_num = 0
        ins_num = del_num = sub_num = 0
        af_code = []

        for i in range(y.size):
            unit_list = [3, 2, 1, 0]
            ins_times = 0
            while ins_times < 1:
                if np.random.uniform(0, 1) <= ins_pr:
                    af_code.append(random.choice(unit_list))
                    ins_num = ins_num + 1
                else:
                    break
            if np.random.uniform(0, 1) <= del_pr:
                del_num += 1
                continue
            else:
                pt_num += 1
                if np.random.uniform(0, 1) <= sub_pr:
                    unit_list.pop(y[i])
                    af_code.append(random.choice(unit_list))
                    sub_num += 1
                else:
                    af_code.append(y[i])

        af_code_array = np.array(af_code)
        return af_code_array, del_num, ins_num, sub_num

    def random_sample(self, sample_times: int):
        """
        Generate random error samples from input sequences.

        Parameters:
        - sample_times: Number of samples to generate.

        Returns:
        - sample_packet: List of sequences with random errors applied.
        """
        pr = self.pr
        strands = self.strands
        sample_packet = []
        dropout_times = 0
        total_ins_num = 0
        total_del_num = 0
        total_sub_num = 0
        ex_error_seq = 0

        for i in range(len(strands)):
            for j in range(sample_times):
                if round(np.random.uniform(0, 1), 5) >= pr:
                    sample_strand = strands[i]
                    error_strand, del_num, ins_num, sub_num = self.random_error(sample_strand)
                    sample_packet.append(error_strand)
                    total_ins_num += ins_num
                    total_del_num += del_num
                    total_sub_num += sub_num

                    if del_num + ins_num + sub_num >= 2:
                        ex_error_seq += 1

                else:
                    dropout_times += 1

        print("Number of sequences with more than two insertions, deletions, or substitutions:", ex_error_seq)
        print("Number of dropped sequences:", dropout_times)
        print("Total number of substitution/deletion/insertion errors:", total_sub_num, total_del_num, total_ins_num)
        random.shuffle(sample_packet)
        return sample_packet


 
       