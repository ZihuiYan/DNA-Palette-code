import numpy as np

class VTCode:
    def __init__(self, n: int, q: int, a = 0, b = 0,
                correct_substitutions = False):
        '''
        Here n is the codeword length and q is the alphabet size.
        a and b are parameters of the code that do not impact the rate in this
        implementation (so can be left at their default values).
        Set correct_substitutions to True for q = 2 if you want ability to correct
        single substitution errors as well.
        '''
        assert q >= 2
        assert n >= 2
        self.n = n
        self.q = q
        self.correct_substitutions = correct_substitutions
        self.k = find_k(self.n, self.q, self.correct_substitutions)
        assert self.k > 0
        self.a = a
        self.b = b
        if self.q == 2:
            if not self.correct_substitutions:
                self.m = self.n + 1
            else:
                self.m = 2*self.n + 1
            assert 0 <= self.a < self.m
            self._generate_systematic_positions_binary()
        else:
            self.m = self.n
            self.t = np.ceil(np.log2(n)).astype(np.int64)
            assert 0 <= self.a < self.m
            assert 0 <= self.b < self.q
            self._generate_tables()
            


    def decode(self, y):
        '''
        input  y: list or 1d np array with the noisy codeword
        return x: decoded message bits as a 1d numpy array with dtype int64 or
                  None if decoding fails
        '''
        y = np.array(y, dtype=np.int64)
        assert y.ndim == 1
        n_y = y.size
        if (n_y < self.n - 1) or (n_y > self.n + 1):
            return None
        if (np.max(y) > self.q-1) or (np.min(y) < 0):
            print("Value in y out of range 0...q-1")
            raise RuntimeError
        if self.q == 2:
            if n_y != self.n:
                y = _correct_binary_indel(self.n, self.m, self.a, y)
            else:
                if self.correct_substitutions and not self._is_codeword(y):
                    y = _correct_binary_substitution(self.n, self.m, self.a, y)
        else:
            if n_y != self.n:
                y = _correct_q_ary_indel(self.n, self.m, self.a, self.b, self.q, y)
        
        return self._decode_codeword(y)


    def encode(self, x):
        '''
        input  x: list or 1d np array with the message bits (length k)
        return y: encoded codeword as a 1d numpy array with dtype int64 (length n)
        '''
        x = np.array(x, dtype = np.int64)
        assert x.ndim == 1
        #assert x.size == self.k
        if (np.max(x) > 1) or (np.min(x) < 0):
            print("Value in x out of range {0, 1}")
            raise RuntimeError
        if self.q == 2:
            return self._encode_binary(x)
        else:
            return self._encode_q_ary(x)


    def _decode_codeword(self, y):
        '''
        decode a codeword (if not a codeword, returns None)
        '''
        if not self._is_codeword(y):
            return None
        if self.q == 2:
            return self._decode_codeword_binary(y)
        else:
            return self._decode_codeword_q_ary(y)

    def _decode_codeword_binary(self, y):
        '''
        decoding helper for binary case (assume it's a valid codeword)
        '''
        # just return values at the systematic positions
        return y[self.systematic_positions-1]



    def _encode_binary(self, x):
        '''
        encoding helper for binary case
        '''
        y = np.zeros(self.n, dtype = np.int64)
        # first set systematic positions
        y[self.systematic_positions-1] = x
        # now set the rest positions based on syndrome
        syndrome = _compute_syndrome_binary(self.m, self.a, y)
        if syndrome != 0:
            for pos in reversed(self.parity_positions):
                if syndrome >= pos:
                    y[pos-1] = 1
                    syndrome -= pos
                    if syndrome == 0:
                        break
        assert self._is_codeword(y)
        return y


    def _is_codeword(self, y):
        '''
        return True if y is a codeword
        '''
        if y is None or y.size != self.n:
            return False
        if self.q == 2:
            return (_compute_syndrome_binary(self.m, self.a, y) == 0)
        else:
            return (_compute_syndrome_q_ary(self.m, self.a, self.b, self.q, y) == (0,0))

    def _generate_systematic_positions_binary(self):
        # generate positions of systematic and parity bits (1 indexed)
        t = np.ceil(np.log2(self.n+1)).astype(np.int64)
        # put powers of two in the parity positions
        self.parity_positions = np.zeros(self.n-self.k, dtype=np.int64)
        for i in range(t):
            self.parity_positions[i] = np.power(2,i)
        if self.correct_substitutions:
            assert self.parity_positions.size == t + 1
            # one extra parity bit in this case
            # depending on if last position in codeword is already filled,
            # set it or the previous position as a parity_position
            if self.parity_positions[t-1] == self.n:
                self.parity_positions[t-1] = self.n - 1
                self.parity_positions[t] = self.n
            else:
                self.parity_positions[t] = self.n
        self.systematic_positions =  np.setdiff1d(np.arange(1,self.n+1), self.parity_positions)
        return



def find_smallest_n(k: int, q : int, correct_substitutions = False):
    '''
    Returns smallest n for a code with given k and q.
    Here k is the message length in bits, n is the codeword length and q is the
    alphabet size.
    Set correct_substitutions to True for q = 2 if you want ability to correct
    single substitution errors as well.
    '''
    assert q >= 2
    assert k >= 1
    if q != 2 and correct_substitutions == True:
        print("correct_substitutions can be True only for q = 2")
        raise RuntimeError
    # set the starting n and then increase till you get the minimum
    if q == 2:
        if not correct_substitutions:
            n = k + np.ceil(np.log2(k+1)).astype(np.int64)
        else:
            n = k + np.ceil(np.log2(2*k+1)).astype(np.int64)
    else:
        n = int(k/np.ceil(np.log2(q)).astype(np.int64))
    while True:
        if find_k(n, q, correct_substitutions) >= k:
            break
        n += 1
    return n

def find_k(n: int, q: int, correct_substitutions = False):
    '''
    Returns k for a code with given n and q.
    Here k is the message length in bits, n is the codeword length and q is the
    alphabet size.
    Set correct_substitutions to True for q = 2 if you want ability to correct
    single substitution errors as well.
    '''
    if q != 2 and correct_substitutions:
        print("correct_substitutions can be True only for q = 2")
        raise RuntimeError
    if q == 2:
        if not correct_substitutions:
            return n - np.ceil(np.log2(n+1)).astype(np.int64)
        else:
            return n - np.ceil(np.log2(2*n+1)).astype(np.int64)
    else:
        t = np.ceil(np.log2(n)).astype(np.int64)
        if q == 3:
            if n < 7:
                return 0
            if _power_of_two(n-1):
                # in this case we can't store data in 2**(t-1)+1
                return np.floor((n-3*t+3)*np.log2(q)).astype(np.int64) + 2*(t-4) + 1
            else:
                return np.floor((n-3*t+3)*np.log2(q)).astype(np.int64) + 2*(t-3)
        else:
            if n < 6:
                return 0
            if _power_of_two(n-1):
                # in this case we can't store data in 2**(t-1)+1
                return np.floor(np.max([(n-3*t+3),0])*np.log2(q)).astype(np.int64) + \
                    np.floor(2*np.log2(q-1)).astype(np.int64)*np.max([(t-4),0]) + \
                    2*np.floor(np.log2(q-1)).astype(np.int64)
            else:
                return np.floor(np.max([(n-3*t+3),0])*np.log2(q)).astype(np.int64) + \
                    np.floor(2*np.log2(q-1)).astype(np.int64)*np.max([(t-3),0]) + \
                    np.floor(np.log2(q-1)).astype(np.int64)


# internal functions
def _correct_binary_indel(n: int, m: int, a: int, y):
    '''
    correct single insertion deletion error in the binary case.
    Used also in the q-ary alphabet case as a subroutine.
    Input: n (codeword length), m (modulus), a (syndrome),
           y (noisy codeword, np array)
    Output: corrected codeword
    '''
    s = _compute_syndrome_binary(m, a, y)
    w = np.sum(y)
    y_decoded = np.zeros(n, dtype=np.int64)
    if y.size == n-1:
        # deletion
        if s == 0:
            # last entry 0 was deleted
            y_decoded[:-1] = y
        elif s <= w:
            # 0 deleted and s = number of 1s to right
            num_ones_seen = 0
            for i in reversed(range(n-1)):
                if y[i] == 1:
                    num_ones_seen += 1
                    if num_ones_seen == s:
                        y_decoded[:i] = y[:i]
                        y_decoded[i+1:] = y[i:]
                        break
        else:
            # 1 deleted and s-w-1 = number of 0s to left
            num_zeros_seen = 0
            if s-w-1 == 0:
                y_decoded[0] = 1
                y_decoded[1:] = y
            else:
                success = False
                for i in range(n-1):
                    if y[i] == 0:
                        num_zeros_seen += 1
                        if num_zeros_seen == s-w-1:
                            y_decoded[:i+1] = y[:i+1]
                            y_decoded[i+1] = 1
                            y_decoded[i+2:] = y[i+1:]
                            success = True
                            break
                if not success:
                    y_decoded = None
    else:
        # insertion
        if s == m-n-1 or s == 0:
            # last entry inserted
            y_decoded = y[:-1]
        elif s == m-w:
            # remove first entry
            y_decoded = y[1:]
        elif s > m-w:
            # 0 was inserted, m-s 1's to the right of this zero
            num_ones_seen = 0
            success = False
            for i in reversed(range(2,n+1)):
                if y[i] == 1:
                    num_ones_seen += 1
                    if num_ones_seen == m-s:
                        if y[i-1] == 0:
                            y_decoded[:i-1] = y[:i-1]
                            y_decoded[i-1:] = y[i:]
                            success = True
                        else:
                            pass
                        break
            if not success:
                y_decoded = None
        else:
            # 1 was inserted, m-w-s 0's to the left of this 1
            num_zeros_seen = 0
            success = False
            for i in range(n-1):
                if y[i] == 0:
                    num_zeros_seen += 1
                    if num_zeros_seen == m-w-s:
                        if y[i+1] == 1:
                            y_decoded[:i+1] = y[:i+1]
                            y_decoded[i+1:] = y[i+2:]
                            success = True
                        else:
                            pass
                        break
            if not success:
                y_decoded = None
    return y_decoded

def _correct_binary_substitution(n: int, m: int, a: int, y):
    '''
    correct single substitution error in the binary case.
    Input: n (codeword length), m (modulus), a (syndrome),
           y (noisy codeword, np array)
    Output: corrected codeword
    '''
    assert m == 2*n+1
    s = _compute_syndrome_binary(m, a, y)
    y_decoded = np.array(y)
    if s == 0:
        # no error, nothing to do
        pass
    elif s < n+1:
        # 1 flipped to 0 at s
        y_decoded[s-1] = 1
    else:
        # 0 flipped to 1 at 2n+1-s
        y_decoded[2*n+1-s-1] = 0
    return y_decoded

def _compute_syndrome_binary(m: int, a: int, y):
    '''
    compute the syndrome in the binary case (a - sum(i*y_i) mod m)
    '''
    n_y = y.size
    return np.mod(a - np.sum((1+np.arange(n_y))*y),m)


def _convert_y_to_alpha(y):
    '''
    convert q-ary y of length n to binary length n-1 alpha, alpha_i = 1 iff y_i >= y_{i-1}
    '''
    return (y[1:] >= y[:-1]).astype(np.int64)


def _convert_base(in_array, in_base, out_base, out_size = None):
    # convert in_array represented in in_base to array in out_base (both MSB first)
    # pad to out_size
    return _number_to_q_ary_array(_q_ary_array_to_number(in_array, in_base), out_base, out_size)

def _power_of_two(num):
    '''
    return True if num is a power of 2
    '''
    return np.ceil(np.log2(num)).astype(np.int64) == np.floor(np.log2(num)).astype(np.int64)
