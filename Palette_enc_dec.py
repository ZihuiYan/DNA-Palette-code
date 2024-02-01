import numpy as np
import math
import random
from vt_enc import vt_encode, vt_decode
from collections import Counter,defaultdict
import copy


def str_to_int(or_str):
    # Convert the string sequence to an integer in the corresponding base
    times_int = 0
    for j in range(len(or_str)):
        times_int += int(or_str[j]) * 2**(len(or_str) - j - 1)
    return times_int

def quanternary_read(binary_string):
    # Convert binary string to quaternary string
    quan_str = ''.join(str(int(binary_string[i:i + 2], 2)) for i in range(0, len(binary_string), 2))
    return quan_str

def quan_to_int(quan_str):
    # Convert quaternary string to integer
    int_value = sum(int(digit) * 4**(len(quan_str) - i - 1) for i, digit in enumerate(quan_str))
    return int_value

def quan_to_binary(quan_str):
    # Convert quaternary string to binary string
    binary_str = ''.join(bin(int(digit, 4))[2:].zfill(2) for digit in quan_str)
    return binary_str

def vote(codeword_list):
    # Majority vote, determine the most common element
    result_list = [Counter(column).most_common(1)[0][0] for column in zip(*codeword_list)]
    return result_list

def xor_encrypt(input_string, key):
    # Encrypt using XOR with the given key
    encrypted = ''.join(str((int(input_char) + int(key_char)) % 4) for input_char, key_char in zip(input_string, key))
    return encrypted

def xor_decrypt(input_string, key):
    # Decrypt using XOR with the given key
    decrypted = ''.join(str((int(input_char) - int(key_char)) % 4) for input_char, key_char in zip(input_string, key))
    return decrypted


def check_function(or_msg, dec_txt):
    """
    Check the number of differences between two strings.

    Parameters:
    - or_msg: Original message.
    - dec_txt: Decoded text.

    Returns:
    - The number of differences.
    """
    err_times = sum(1 for a, b in zip(or_msg, dec_txt) if a != b)
    if len(or_msg) != len(dec_txt):
        print("Mismatch in lengths!", len(or_msg), len(dec_txt))
    return err_times


def Palette_enc(binary_string, length_ary, id_num, toatl_file_num, key_str):
    """
    Encode the given binary string using a specified method.

    Parameters:
    - binary_string: Original binary string.
    - length_ary: Radix.
    - id_num: ID of the current file.
    - toatl_file_num: Total number of files.
    - key_str: Seed to determine the random key string.

    Returns:
    - Codeword table, minus sequence number, length of the original binary string, code word length, and number of VT information bits.
    """
    enc_table = []
    minus_num = 0
    b_len = math.ceil(math.log2(len(binary_string) / length_ary))  # Number of label positions
    if b_len % 2 != 0:
        b_len += 1

    id_len = math.ceil(math.log2(toatl_file_num))
    if id_len % 2 != 0:
        id_len += 1
    id_index = bin(id_num)[2:].rjust(id_len, '0')

    quan_id_index = quanternary_read(id_index)  # Current file index

    if length_ary == 272: # Key word for DICOM
        key = '02131213030120302130320121210321303132020231031021312031212021310123203010213120203012230303120121311032110320123020203131023031230121320301202132'
    else:
        random.seed(key_str)
        key_length = (id_len + b_len + length_ary) // 2
        key = ''.join(random.choice('0123') for _ in range(key_length)) # Key word for others

    ab_len = (length_ary // 2) // (b_len // 2)
    result_str = False
    for i in range((len(binary_string) // length_ary)):
        a_str = binary_string[length_ary * i:length_ary * (i + 1)]

        if all(c == '0' for c in a_str):
            minus_num += 1
        else:
            a_quan_str = quanternary_read(a_str)
            b_str = bin(i)[2:].rjust(b_len, '0')
            b_quan_str = quanternary_read(b_str)

            ab_quan_str = ''
            for t in range(b_len // 2):
                ab_quan_str += a_quan_str[ab_len * t:ab_len * (t + 1)]
                ab_quan_str += b_quan_str[t]
            ab_quan_str += a_quan_str[ab_len * (b_len // 2):]

            word_ab = ''
            for t in range(len(quan_id_index)):
                word_ab += quan_id_index[t]
                word_ab += ab_quan_str[
                    (len(ab_quan_str) // len(quan_id_index)) * t:(len(ab_quan_str) // len(quan_id_index)) * (t + 1)]
            word_ab += ab_quan_str[(len(ab_quan_str) // len(quan_id_index)) * len(quan_id_index):]

            word_key = xor_encrypt(word_ab, key)
            code_word, code_word_len0, k0 = vt_encode(word_key)
            enc_table.append(code_word)
            if not result_str:
                a_len1 = len(a_quan_str)
                b_len1 = len(b_quan_str)
                code_word_len1 = code_word_len0
                k1 = k0
                result_str = True

    if len(binary_string) % length_ary != 0:
        a_str = binary_string[length_ary * (len(binary_string) // length_ary):]
        if all(c == '0' for c in a_str):
            minus_num += 1
        else:
            a_quan_str = quanternary_read(a_str).ljust(a_len1, '0')
            b_str = bin((len(binary_string) // length_ary))[2:].rjust(b_len, '0')
            b_quan_str = quanternary_read(b_str)

            ab_quan_str = ''
            for i in range(len(a_quan_str)):
                if i < len(b_quan_str):
                    ab_quan_str += b_quan_str[i]
                    ab_quan_str += a_quan_str[i]
                else:
                    ab_quan_str += a_quan_str[i]
            word_ab = quan_id_index + ab_quan_str
            word_key = xor_encrypt(word_ab, key)
            code_word, code_word_len0, k0 = vt_encode(word_key)
            enc_table.append(code_word)
    return enc_table, minus_num, len(binary_string), code_word_len1, k1


def Palette_dec(txt_table, length_ary, enc_binary_len_dic, code_word_len, vt_k, toatl_file_num, key_str):
    """
    Decode the given table of text.

    Parameters:
    - txt_table: Table where each entry is a string representing a sequence.
    - length_ary: Radix.
    - enc_binary_len_dic: Dictionary containing the length of the original binary sequence for each file.
    - code_word_len: Length of the code word.
    - vt_k: VT information bits.
    - toatl_file_num: Total number of files.
    - key_str: Key string.

    Returns:
    - Dictionary containing the decoded binary sequence for each file.
    """

    id_len = math.ceil(math.log2(toatl_file_num))
    if id_len % 2 != 0:
        id_len += 1

    error_seq_num = 0

    if length_ary == 272:
        key = '02131213030120302130320121210321303132020231031021312031212021310123203010213120203012230303120121311032110320123020203131023031230121320301202132'
    else:
        random.seed(key_str)
        key = ''.join(random.choice('0123') for _ in range(vt_k))

    file_dic = {}
    dec_dic = {}

    for i in range(len(txt_table)):
        code_word0 = vt_decode(txt_table[i], code_word_len, vt_k)
        ab_quan_len = (vt_k - id_len // 2) // (id_len // 2)

        if code_word0 is not None:
            code_word0 = xor_decrypt(code_word0, key)
            id_index = ''
            ab_quan_str = ''

            for _ in range(id_len // 2):
                id_index += code_word0[0]
                ab_quan_str += code_word0[1:1 + ab_quan_len]
                code_word0 = code_word0[1 + ab_quan_len:]
            ab_quan_str += code_word0
            id_num = quan_to_int(id_index)

            if id_num <= toatl_file_num - 1:
                file_dic.setdefault(id_num, []).append(ab_quan_str)

        else:
            error_seq_num += 1

    for id_num in file_dic.keys():
        enc_binary_len = enc_binary_len_dic[id_num]
        codeword_dic = defaultdict(list)
        dec_txt = [0] * enc_binary_len
        b_len = math.ceil(math.log2(enc_binary_len / length_ary))
        if b_len % 2 != 0:
            b_len += 1
        ab_len = (length_ary // 2) // (b_len // 2)

        for code_word in file_dic[id_num]:
            b_quan_str = ''
            a_quan_str = ''
            code_word_ab = copy.copy(code_word)

            for _ in range(b_len // 2):
                a_quan_str += code_word_ab[0:ab_len]
                b_quan_str += code_word_ab[ab_len]
                code_word_ab = code_word_ab[ab_len + 1:]
            a_quan_str += code_word_ab

            b_index = quan_to_int(b_quan_str)
            a_binary_str = quan_to_binary(a_quan_str)
            codeword_dic[b_index].append(a_binary_str)

        for b_index, a_str_list in codeword_dic.items():
            a_str = vote(a_str_list)

            if length_ary * (b_index + 1) <= enc_binary_len:
                dec_txt[length_ary * b_index:length_ary * (b_index + 1)] = a_str
            else:
                dec_txt[length_ary * b_index:] = a_str[:enc_binary_len - (length_ary * b_index)]
        dec_bianry = "".join('%s' % id for id in dec_txt)
        dec_dic[id_num] = dec_bianry

    #print("Number of failed VT decodings:", error_seq_num)
    return dec_dic

    