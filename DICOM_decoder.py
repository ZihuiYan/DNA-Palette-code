from DNA_Ladder_code import differential_dec, readFile, read_file, rs_decode
from Palette_enc_dec import check_function, Palette_dec


# File paths
parent_dir = "./DICOM_files"
file_path_list = readFile(parent_dir)
print(file_path_list)


# Input file paths and output dictionary
input_file_paths = ["sequence reads/sequences_0.0185.txt"]
output_dict = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 0}
total_codeword_table = []

# Process input files
for input_file_path in input_file_paths:
    with open(input_file_path, 'r') as f:
        for line in f:
            code_list0 = [output_dict[base] for base in line.strip()]
            total_codeword_table.append(code_list0)

print("Sequencing Reads number:", len(total_codeword_table))
print("Average Coverage: {:03f} ".format(len(total_codeword_table) / 255248))

# Encoding information
enc_binary_len_dic = {0: 2745568, 1: 2745568, 2: 2745568, 3: 2745568, 4: 2745568, 5: 2745840, 6: 2745568, 7: 2745568,
                      8: 2745568, 9: 2745568, 10: 2745568, 11: 2745568, 12: 2745568, 13: 2745568, 14: 2745568, 15: 2745568,
                      16: 2745568, 17: 2745568, 18: 2745568, 19: 2745568, 20: 2745568, 21: 2409104, 22: 2409104, 23: 2409104,
                      24: 2409104, 25: 2409104, 26: 2409104, 27: 2409104, 28: 2409104, 29: 2409104, 30: 2409104, 31: 2409104,
                      32: 2409104, 33: 2409104, 34: 2409104, 35: 2409104, 36: 2409104, 37: 2409104, 38: 2409104, 39: 2409104,
                      40: 2409104, 41: 2409104}
file_len_dic = {0: 299670, 1: 299674, 2: 299672, 3: 299676, 4: 299674, 5: 299678, 6: 299670, 7: 299674, 8: 299670,
                9: 299672, 10: 299672, 11: 299670, 12: 299670, 13: 299672, 14: 299672, 15: 299674, 16: 299674, 17: 299672,
                18: 299672, 19: 299672, 20: 299674, 21: 263042, 22: 263042, 23: 263042, 24: 263042, 25: 263042, 26: 263044,
                27: 263044, 28: 263044, 29: 263044, 30: 263044, 31: 263046, 32: 263044, 33: 263044, 34: 263044, 35: 263044,
                36: 263042, 37: 263040, 38: 263040, 39: 263040, 40: 263040, 41: 263040}
str_len = 136

# Decoding stage
dec_dic = Palette_dec(total_codeword_table, str_len * 2, enc_binary_len_dic, 155, 146, 42, 0)
rs_dec_dic = rs_decode(dec_dic, file_len_dic, str_len // 4)

# Verification stage
total_byte_num = 0
total_file_byte = 0
correct_file = 0
total_error_num = 0
for id_num in range(len(file_path_list)):
    rs_dec_file_path = "./rs_dec_files/DICOM{:03}".format(id_num) + '.txt'
    with open(rs_dec_file_path, 'wb') as f:
        f.write(rs_dec_dic[id_num])

total_byte_num = 0
total_file_byte = 0
correct_file = 0
total_error_num = 0
Ladder_file_path_list = readFile("./rs_dec_files")

re_file_path = "./dec_files"
differential_dec(Ladder_file_path_list, re_file_path)
re_file_list = readFile(re_file_path)

for id_num in range(len(file_path_list)):
    dec_bytes = read_file(re_file_list[id_num])
    file_bytes = read_file(file_path_list[id_num])
    error_byte = check_function(file_bytes, dec_bytes)
    total_file_byte += len(file_bytes)
    total_error_num += error_byte
    byte_error_rate = error_byte / len(file_bytes)
    print("Byte error rate for {}-th file: {:.3f}".format(id_num, byte_error_rate))
    if error_byte == 0:
        correct_file += 1

print("Byte error rate: {:03f}".format(total_error_num / total_file_byte))
print("Correctly recovered file num: ", correct_file)

