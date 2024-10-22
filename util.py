from config import *
import numpy as np

def parse_seqs_file(fname):
    with open(fname, 'r') as seq_file:
        allSeqs = []
        numSeq = 0
        for tmp_str in seq_file:
            tmp_str = tmp_str.strip()
            ht_tmp_seq = ['*'] + list(tmp_str) + ['#']
            allSeqs.append(ht_tmp_seq)
            numSeq += 1
    return allSeqs, numSeq

def get_init_model_length(lenSeqs):
    return max(lenSeqs)

def sequence_dump(allSeqs):
    for seq in allSeqs:
        print(''.join(seq))

def identity_4D(y):
    I1 = np.eye(y.shape[0])
    I2 = np.eye(y.shape[1])
    I3 = np.eye(y.shape[2])
    I4 = np.eye(y.shape[3])
    big_eye = np.zeros((sum(y.shape),sum(y.shape)))
    big_eye[:y.shape[0],:y.shape[0]] = I1
    big_eye[y.shape[0]:y.shape[0]+y.shape[1], y.shape[0]:y.shape[0]+y.shape[1]] = I2
    big_eye[y.shape[0]+y.shape[1]:y.shape[0]+y.shape[1]+y.shape[2], y.shape[0]+y.shape[1]:y.shape[0]+y.shape[1]+y.shape[2]] = I3
    big_eye[y.shape[0]+y.shape[1]+y.shape[2]:, y.shape[0]+y.shape[1]+y.shape[2]:] = I4
    return big_eye

def print_cube(cube):
    T1 = len(cube)
    T2 = len(cube[0])
    T3 = len(cube[0][0])
    for i in range(T1):
        for j in range(T2):
            for k in range(T3):
                print(cube[i][j][k].toString())

def Frobenius_prod(A, B):
    p = A * B
    return np.sum(p, axis=None)

def print_C(C):
    T0 = len(C)
    for n in range(T0):
        T1, T2, T3, T4 = C[n].shape
        for i in range(T1):
            for j in range(T2):
                for k in range(T3):
                    for m in range(T4):
                        print(f"C[{n}][{i}][{j}][{k}][{m}] = {C[n][i][j][k][m]}")

def move2T3idx(m):
    mapping = {
        0: -1,
        DELETION_A: 0,
        MATCH_A: 0,
        DELETION_T: 1,
        MATCH_T: 1,
        DELETION_C: 2,
        MATCH_C: 2,
        DELETION_G: 3,
        MATCH_G: 3,
        DELETION_START: 4,
        MATCH_START: 4,
        DELETION_END: 5,
        MATCH_END: 5
    }
    return mapping.get(m, -2)

def dna2T3idx(dna):
    m = {'A': 0, 'T': 1, 'C': 2, 'G': 3, '*': 4, '#': 5, GAP_NOTATION: -1}
    if dna in m:
        return m[dna]
    else:
        print("error", dna)
        exit()

def T3idx2dna(idx):
    m = ['A', 'T', 'C', 'G', '*', '#']
    if idx < len(m) and idx >= 0:
        return m[idx]
    else:
        print("error")
        exit()
        return 

action2str = {
        INSERTION: "Insertion",
        DELETION_A: "Deletion_A",
        DELETION_T: "Deletion_T",
        DELETION_C: "Deletion_C",
        DELETION_G: "Deletion_G",
        DELETION_START: "DELETION_START",
        DELETION_END: "DELETION_END",
        MATCH_A: "Match_A",
        MATCH_T: "Match_T",
        MATCH_C: "Match_C",
        MATCH_G: "Match_G",
        MATCH_START: "MATCH_START",
        MATCH_END: "MATCH_END",
        UNDEFINED: "Undefined"
    }

class Cell:
    def __init__(self, dim):
        self.dim = dim
        self.score = 0
        self.action = UNDEFINED
        self.location = [-1 for _ in range(self.dim)]
        self.acidA = '?'
        self.acidB = '?'
        self.ans_idx = -1

    def toString(self):
        s = []
        s.append("[(")
        for i in range(self.dim):
            s.append(self.location[i])
            if (i < self.dim - 1):
                s.append(",")

        s.append("), ")
        s.append(action2str.get(self.action))
        s.append(", ")
        s.append(self.acidA)
        s.append(", ")
        s.append(self.acidB)
        s.append(", ")
        s.append(f'{self.score:g}')
        s.append(") ")
        return ''.join([str(x) for x in s])
    