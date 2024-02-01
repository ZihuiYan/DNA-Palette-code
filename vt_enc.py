#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 20:39:00 2023

@author: serena-mo
"""

from vt import VTCode,find_smallest_n
import numpy as np
import ast
import re


def deleaved_str(msg):
    msg_0 = np.zeros(len(msg))
    msg_1 = np.zeros(len(msg))
    for i in range(len(msg)):
        msg_0[i] = (msg[i] % 2)
        msg_1[i] = (msg[i] // 2)
    return msg_0,msg_1

def leaved_str(msg_0,msg_1):
    msg = np.zeros(len(msg_0))
    for i in range(len(msg_0)):
        msg[i] = msg_0[i] + msg_1[i] * 2
    return msg

def vt_encode(msg):
    msg_list = [int(char) for char in msg]
    #print(msg_list)
    msg_array = np.array(msg_list)
    #print(msg_array)
    n = find_smallest_n(len(msg_array), 2,correct_substitutions = True)
    msg_0,msg_1 = deleaved_str(msg_array)
    enc_msg_0 = VTCode(n, 2, 0, 0, correct_substitutions = True).encode(msg_0)
    enc_msg_1 = VTCode(n, 2, 0, 0, correct_substitutions = True).encode(msg_1)
    enc_msg = leaved_str(enc_msg_0, enc_msg_1)
    enc_msg = np.array(enc_msg,dtype=np.int64)
    return enc_msg, n, len(msg_array)

def vt_decode(channel_output, n, k):
    code = VTCode(n, 2, 0, 0, correct_substitutions = True)
    enc_0,enc_1 = deleaved_str(channel_output)
    
    dec_0 = code.decode(enc_0)    
    dec_1 = code.decode(enc_1)

    if dec_1 is None:
        return None
    elif dec_0 is None:
        return None
    else:   
        dec_word = leaved_str(dec_0,dec_1).astype(int)
        dec_str = ''.join(str(i) for i in dec_word)
        #print(type(dec_str))
        return dec_str

if __name__ == '__main__':
    msg = '20032101212111111031111203230320023303003011023303030300022202130000011'
    enc, n, k = vt_encode(msg)
    print(enc)
    dec = vt_decode(enc, n, k)
    print(dec,len(dec))