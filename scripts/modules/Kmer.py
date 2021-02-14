import numpy
import math


class BaseConverter:
    def __init__(self):
        self.base_to_index = {
            "A":0,
            "C":1,
            "G":2,
            "T":3
        }

        self.index_to_base = ["A","C","G","T"]

        self.byte_to_index = numpy.array([
            4,4,4,4,4,4,4,4,4,4,      # 0
            4,4,4,4,4,4,4,4,4,4,      # 10
            4,4,4,4,4,4,4,4,4,4,      # 20
            4,4,4,4,4,4,4,4,4,4,      # 30
            4,4,4,4,4,4,4,4,4,4,      # 40
            4,4,4,4,4,4,4,4,4,4,      # 50
            4,4,4,4,4,0,4,1,4,4,      # 60  A = 65, C = 67
            4,2,4,4,4,4,4,4,4,4,      # 70  G = 71
            4,4,4,4,3,4,4,4,4,4,      # 80  T = 84
            4,4,4,4,4,4,4,4,4,4,      # 90
            4,4,4,4,4,4,4,4,4,4,      # 100
            4,4,4,4,4,4,4,4,4,4,      # 110
            4,4,4,4,4,4,4,4,4,4,      # 120
        ], dtype=numpy.uint8)


class Kmer:
    def __init__(self, k):
        self.k = k
        self.number = 0
        self.base_converter = BaseConverter()

    def __str__(self):
        n = self.number
        s = ""

        mask = 2**(2*self.k-1) + 2**(2*self.k - 2)

        for i in range(self.k):
            s += self.base_converter.index_to_base[(n & mask) >> 2*(self.k - i - 1)]
            mask >>= 2

        return s


class KmerTable:
    def __init__(self, k, density):
        self.k = k
        self.n_markers = 0
        self.density = density
        self.threshold = int(round(float(101)*float(density)))
        self.max_value = 4**self.k

        # Mapping from bitwise kmer values (binary 2 bit nucleotides) to table index
        self.table = numpy.full(self.max_value, fill_value=-1, dtype=numpy.int64)
        self.flag_kmers()

    def flag_kmers(self):
        self.n_markers = 0

        for i in range(self.max_value):
            if i % 101 < self.threshold:
                self.table[i] = self.n_markers
                self.n_markers += 1

    def at(self, kmer_number):
        return self.table[kmer_number]


class KmerIterator:
    def __init__(self, sequence, k):
        self.sequence = sequence
        self.k = k
        self.i = 0
        self.kmer = Kmer(k)
        self.base_converter = BaseConverter()

        self.initialize_kmer()

        if self.k > len(self.sequence):
            exit("ERROR: k (%d) is longer than sequence" + str(self.sequence))

    def next(self):
        self.kmer.number <<= 2

        if self.i >= len(self.sequence):
            return False

        if type(self.sequence) == str:
            base_index = self.base_converter.byte_to_index[ord(self.sequence[self.i])]
        else:
            base_index = self.base_converter.byte_to_index[self.sequence[self.i]]

        if base_index < 4:
            self.kmer.number |= base_index
            self.kmer.number &= 4**self.k - 1
        else:
            exit("ERROR: base not found in 2bit translation table (not ACGT)")

        self.i += 1

        return True

    def initialize_kmer(self):
        for i in range(self.k - 1):
            self.next()

        return


class MarkerIterator:
    def __init__(self, sequence, k, density):
        self.kmer_iterator = KmerIterator(sequence, k=k)
        self.density = density
        self.k = k
        self.threshold = int(round(float(101)*float(density)))

    def next(self):
        success = False

        while self.kmer_iterator.next():
            if self.kmer_iterator.kmer.number % 101 < self.threshold:
                success = True
                break

        return success

    def kmer_number(self):
        return self.kmer_iterator.kmer.number

    def marker_number(self):
        return int((round(float(self.kmer_iterator.kmer.number)/101))*(self.density*100) + (self.kmer_iterator.kmer.number%101))

    def kmer_string(self):
        return str(self.kmer_iterator.kmer)


def test():
    sequence = "AAAACCCCGGGGTTTTAACACTATGAGACATAGGGTGCTCTGCGGCTGGCTGCGGGTCTTTTCCCC"
    k = 10
    density = 0.2

    kmer_table = KmerTable(k, density)

    print("n_markers", kmer_table.n_markers)

    kmer_iterator = KmerIterator(sequence=sequence, k=k)

    while kmer_iterator.next():
        print(kmer_iterator.kmer, kmer_iterator.kmer.number, "{0:08b}".format(kmer_iterator.kmer.number))

    marker_iterator = MarkerIterator(sequence, k=k, density=density)

    while marker_iterator.next():
        print(marker_iterator.kmer_string(), marker_iterator.kmer_number(), kmer_table.at(marker_iterator.kmer_number()), marker_iterator.marker_number())


if __name__ == "__main__":
    test()
