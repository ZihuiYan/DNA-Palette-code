from DNA_Ladder_code import differential_enc, readFile, read_file, rs_encode
from Palette_enc_dec import Palette_enc

# Differential encoding
parent_dir = "./DICOM_files"
file_path_list = readFile(parent_dir)
preprocess_file = "./data_preprocess_files"
differential_enc(file_path_list, preprocess_file)

file_path_list = readFile(preprocess_file)
print(file_path_list)

# Initialize evaluation metrics
total_error_num = 0
total_byte_num = 0
total_seq_num = 0
total_nt_num = 0
total_minus_seq = 0
str_len = 136  # Length of the a sequence
total_codeword_table = []
file_len_dic = {}
enc_binary_len_dic = {}
code_word_len_dic = {}
k_dic = {}
key_word = 1

# Encoding stage
for i in range(0, len(file_path_list)):
    print('*****Encoding File {:d}******'.format(i))
    # Read information, output byte representation of the original information with XOR reference
    msg_or = read_file(file_path_list[i])
    print('Original information sequence length:', len(msg_or))
    total_byte_num += len(msg_or)
    print('-----Encoding-----')
    # Perform RS encoding, output byte stream with RS check bits added
    rs_enc, file_txt_len = rs_encode(msg_or, str_len // 4)
    # Perform Palette encoding, output codewords
    codeword_table = Palette_enc(rs_enc, str_len * 2, i, len(file_path_list), key_word)
    total_minus_seq += codeword_table[1]
    total_seq_num += len(codeword_table[0])
    total_nt_num += len(codeword_table[0]) * len(codeword_table[0][0])
    total_codeword_table.extend(codeword_table[0])

    # Input the codewords into a new file
    with open('./codeword.txt', 'a') as f:
        output_dict = {0: "A", 1: "T", 2: "G", 3: "C"}
        for t in range(len(codeword_table[0])):
            dna_str = 'CCACGCGTACCGATAGCTTCAG'  # Primer
            for j in range(len(codeword_table[0][t])):
                dna_str += output_dict[codeword_table[0][t][j]]
            dna_str += 'GCAATTGACCCACGCATGTATC'  # Primer
            f.write(dna_str)
            f.write('\n')

print("Length of oligos:", len(total_codeword_table[0]))
print("Total number of encoded bytes:", total_byte_num)
print("Total number of encoded oligos:", total_seq_num)
print("Theoretical number of encoded sequences:", total_seq_num + total_minus_seq)
print("Information density (bits/nt):", total_byte_num * 8 / total_nt_num)
print('------------------------------------------')

'''
# Error simulation stage
error = Error(total_codeword_table, 0.01, 0.00001, 0.00001, 0.00001)
error_seq_table = error.random_sample(1)
'''

