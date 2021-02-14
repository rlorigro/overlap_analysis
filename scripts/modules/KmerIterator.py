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


class KmerTable:
    def __init__(self, k, density):
        self.k = k
        self.density = density
        self.max_value = int(math.ceil(self.density*(4**self.k)))

        # Mapping from bitwise kmer values (binary 2 bit nucleotides) to table index
        self.table = numpy.zeros(self.max_value, dtype=numpy.bool)

    def flag_kmers(self):
        for i in range(self.max_value):
            if int(round(i % 101)) < 10*self.density:
                self.table[i] = 1


class Kmer:
    def __init__(self):
        self.n = 0
        self.base_converter = BaseConverter()

    def __str__(self):


        return


class KmerIterator:
    def __init__(self, sequence, k):
        self.sequence = sequence
        self.k = k
        self.i = 0
        self.kmer = 0
        self.base_converter = BaseConverter()

        self.initialize_kmer()

    def step(self):
        self.kmer <<= 2

        if self.i >= len(self.sequence):
            return False

        base_index = self.base_converter.byte_to_index[ord(self.sequence[self.i])]

        if base_index < 4:
            self.kmer |= base_index
            self.kmer &= 4**self.k - 1
        else:
            exit("ERROR: base not found in 2bit translation table (not ACGT)")

        self.i += 1

        return True

    def initialize_kmer(self):
        for i in range(self.k):
            self.step()

        return


if __name__ == "__main__":
    sequence = "AAAACCCCGGGGTTTT"

    k = 4

    kmer_table = KmerTable(4, 0.1)
    kmer_iterator = KmerIterator(sequence=sequence, k=k)

    print(bin(4**k - 1))

    while(kmer_iterator.step()):
        print(kmer_iterator.kmer, bin(kmer_iterator.kmer))

