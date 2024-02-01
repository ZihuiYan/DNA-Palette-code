import os
import sys
import math
import numpy as np
from reedsolo import RSCodec, ReedSolomonError

import os
import sys


def readFile(parent_path):
    """
    Read a folder and return a list of file names in default order.

    Args:
    - parent_path (str): Folder name.

    Returns:
    - list: List of file names in default order.
    """
    file_path_list = []
    file_list = os.listdir(parent_path)

    for file in file_list:
        if file != '.DS_Store':
            file_path = os.path.join(parent_path, file)
            file_path_list.append(file_path)

    file_path_list.sort()
    return file_path_list


def read_file(file_path: str):
    """
    Read a file and return data in bytearray format.

    Args:
    - file_path (str): File name.

    Returns:
    - bytearray: Data in bytearray format.
    """
    try:
        with open(file_path, "rb") as byte_file:
            byte_string = bytearray(byte_file.read())
        return byte_string
    except OSError:
        print("File not accessible")
        sys.exit()


def dicom_xor_reference(byte_seq0, byte_seq1):  # Input: two byte sequences
    """
    XOR two byte sequences.

    Args:
    - byte_seq0 (bytes): Reference byte sequence.
    - byte_seq1 (bytes): Target byte sequence.

    Returns:
    - bytes: XOR result.
    """
    if len(byte_seq0) < len(byte_seq1):
        add_length = bytearray([0] * (len(byte_seq1) - len(byte_seq0)))
        byte_seq0 = byte_seq0[0:1664] + add_length + byte_seq0[1664:]
    else:
        minus_len = len(byte_seq0) - len(byte_seq1)
        byte_seq0 = byte_seq0[0:1664] + byte_seq0[1664 + minus_len:]

    xor_byte = bytes([a ^ b for a, b in zip(byte_seq0, byte_seq1)])
    return xor_byte


def xor_reference(byte_seq0, byte_seq1):  # Input: two byte sequences
    """
    XOR two byte sequences.

    Args:
    - byte_seq0 (bytes): Reference byte sequence.
    - byte_seq1 (bytes): Target byte sequence.

    Returns:
    - bytes: XOR result.
    """
    if len(byte_seq0) < len(byte_seq1):
        add_length = bytearray([0] * (len(byte_seq1) - len(byte_seq0)))
        byte_seq0 = byte_seq0 + add_length
    else:
        minus_len = len(byte_seq0) - len(byte_seq1)
        byte_seq0 = byte_seq0[minus_len:]

    xor_byte = bytes([a ^ b for a, b in zip(byte_seq0, byte_seq1)])
    return xor_byte


def differential_enc(file_path_list, xor_file_name):
    """
    Perform incremental encoding on files.

    Args:
    - file_path_list (list): List of file names.
    - xor_file_name (str): Output folder for encoded files.
    """
    file_path_list.sort()

    for i in range(0, len(file_path_list)):
        xor_file_path = os.path.join(xor_file_name+"/preprocess_DICOM{:03}".format(i) + '.txt')
        with open(xor_file_path, 'wb') as f:
            if i == 0 or i == 21:
                f.write(read_file(file_path_list[i]))
            else:
                file0 = read_file(file_path_list[i - 1])
                file1 = read_file(file_path_list[i])
                differential_file = dicom_xor_reference(file0, file1)
                f.write(differential_file)


def differential_dec(xor_path_list, re_xor_file_name):
    """
    Restore incrementally encoded files.

    Args:
    - xor_path_list (list): List of file names.
    - re_xor_file_name (str): Output folder for restored files.
    """
    xor_path_list.sort()

    for i in range(0, len(xor_path_list)):
        re_xor_file_path = os.path.join(re_xor_file_name+"/re_DICOM{:03}".format(i) + '.dcm')
        with open(re_xor_file_path, 'wb') as f:
            if i == 0 or i == 21:
                pre_file = read_file(xor_path_list[i])
                f.write(pre_file)
            else:
                file1 = read_file(xor_path_list[i])
                file0 = read_file(re_xor_file_name+"/re_DICOM{:03}".format(i - 1) + '.dcm')
                pre_file = dicom_xor_reference(file0, file1)
                f.write(pre_file)


def rs_encode(file_txt0, a_len):
    """
    RS encode the given file text.

    Parameters:
    - file_txt0: Original file text to be encoded.
    - a_len: Radix.

    Returns:
    - Codeword for DNA Ladder code and original file text length.
    """
    e = 8
    n = 2 ** e - 1
    k = 32
    rsc = RSCodec(k, nsize=n)

    file_txt_len = len(file_txt0)

    file_txt_enc_len = int(math.ceil(file_txt_len / a_len) * a_len)
    if len(file_txt0) % a_len != 0:
        add_bytearray = bytearray(int(a_len - len(file_txt0) % a_len))
        file_txt1 = file_txt0 + add_bytearray
    else:
        file_txt1 = file_txt0

    file_txt_arr = np.array(file_txt1).reshape(file_txt_enc_len // a_len, a_len).T
    enc_word_array = np.zeros(
        (a_len, int(np.ceil((file_txt_enc_len // a_len) / (n - k)) * k + file_txt_enc_len // a_len)), dtype=int)

    for i in range(a_len):
        enc_word_array[i, :] = rsc.encode(file_txt_arr[i])

    enc_word_array = enc_word_array.transpose().flatten()
    enc_word = ''.join(bin(val)[2:].rjust(e, '0') for val in enc_word_array)

    return enc_word, file_txt_len


def rs_decode(dec_dic, file_len_dic, a_len):
    """
    Decode the given dictionary of RS-encoded sequences.

    Parameters:
    - dec_dic: Dictionary containing RS encoded binary sequences for each file.
    - file_len_dic: Dictionary containing the original length of each file.
    - a_len: Radix.

    Returns:
    - Dictionary containing the RS-decoded binary sequence for each file.
    """

    e = 8
    n = 2 ** e - 1
    k = 32
    rsc = RSCodec(k, nsize=n)
    rs_error = 0

    rs_dec_dic = {}

    for id_num in dec_dic.keys():
        rs_binary_str = dec_dic[id_num]
        err_array = bytearray(int(rs_binary_str, 2).to_bytes(len(rs_binary_str) // 8, byteorder='big'))
        dec_txt_arr = np.array(err_array).reshape(len(err_array) // a_len, a_len).T
        dec_word_array = np.zeros((a_len, len(err_array) // a_len - int(np.ceil(len(err_array) / a_len / n) * k)),
                                  np.uint8)

        for t in range(a_len):
            for i in range(len(dec_txt_arr[t]) // n):
                try:
                    dec_word = rsc.decode(dec_txt_arr[t][n * i: n * (i + 1)])[0]
                    dec_word_array[t][(n - k) * i:(n - k) * (i + 1)] = np.array(dec_word)
                except:
                    dec_word_array[t][(n - k) * i:(n - k) * (i + 1)] = dec_txt_arr[t][n * i: n * i + n - k]
                    rs_error += 1

            if len(dec_txt_arr[t]) % n != 0:
                try:
                    dec_word_array[t][(n - k) * (len(dec_txt_arr[t]) // n):] = np.array(
                        rsc.decode(dec_txt_arr[t][n * (len(dec_txt_arr[t]) // n):])[0])
                except:
                    dec_word_array[t][(n - k) * (len(dec_txt_arr[t]) // n):] = dec_txt_arr[t][
                                                                               n * (len(dec_txt_arr[t]) // n): -k]

        dec_word0 = np.copy((dec_word_array.transpose()).flatten())
        dec_word0 = dec_word0[:file_len_dic[id_num]]
        rs_dec_dic[id_num] = dec_word0

    #print('Number of RS decoding errors:', rs_error)
    return rs_dec_dic


if __name__ == '__main__':
    #读取原始文件夹的文件列表并增量编码
    parent_dir = ""
    file_path_list = readFile(parent_dir)
    print(len(file_path_list))
    xor_file_name = ""
    xor_files = differential_enc(file_path_list, xor_file_name)

    #增量解码
    parent_dir = ""
    file_path_list = readFile(parent_dir)
    print(len(file_path_list),file_path_list)
    re_xor_file_name = ""
    re_xor_files = differential_dec(file_path_list, re_xor_file_name)
