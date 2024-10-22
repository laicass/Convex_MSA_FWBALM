import sys, os, time, subprocess
import util
from config import *
from msa_solver import CVX_ADMM_MSA, smith_waterman


def solve(fname):
    allSeqs, numSeq = util.parse_seqs_file(fname)
    lenSeqs = [len(seq) for seq in allSeqs]
    T2 = util.get_init_model_length (lenSeqs) + LENGTH_OFFSET
    util.sequence_dump(allSeqs)
    solver = CVX_ADMM_MSA(allSeqs, lenSeqs, T2)

    print(">>>>>>>>>>>>>>>>>>>>>>>SequenceView<<<<<<<<<<<<<<<<<<<<<<<<")
    recSeq = solver.recSeq
    for n,seq in enumerate(allSeqs):
        print("{:<10}".format("Seq "+str(n)+":"), ''.join(seq))
    print("{:<10}".format("SeqRecov:"), ''.join(recSeq))

    if WRITEFILE:
        fname = os.path.splitext(os.path.basename(fname))[0]
        f = open(fname+".rec", "w")
        f.write(''.join(recSeq[1:-1])+'\n')
        f.close()
    
def local_alignment_cpp(fname):
    fname = os.path.splitext(os.path.basename(fname))[0]
    args = ["MSA_Convex_Local_Pair.exe", fname+".msa", fname+".rec"]
    subprocess.call(args)
    

def local_alignment(allSeqs, numSeq, recSeq):
    #print(">>>>>>>>>>>>>>>>>>>>>>>MatchingView<<<<<<<<<<<<<<<<<<<<<<<<")
    allModelSeqs, allDataSeqs = [], []
    model_seq = recSeq[1:-1]
    for data_seq in allSeqs:
        data_seq = data_seq[1:-1]
        trace = smith_waterman(model_seq, data_seq)
        model_seq = [t.acidA for t in trace]
        data_seq = [t.acidB for t in trace]
        allModelSeqs.append(model_seq)
        allDataSeqs.append(data_seq)
        #print(''.join(model_seq))
        #print(''.join(data_seq))

    print(">>>>>>>>>>>>>>>>>>>>>ClustalOmegaView<<<<<<<<<<<<<<<<<<<<<<")
    allCOSeqs = [[] for _ in range(numSeq)]
    pos = [0 for _ in range(numSeq)]
    while True:
        insertion_ids = set()
        for i in range(numSeq):
            if pos[i] >= len(allModelSeqs[i]):
                continue
            if allModelSeqs[i][pos[i]] == '-':
                insertion_ids.add(i)

        if insertion_ids: # Insertion exists
            for i in range(numSeq):
                if i not in insertion_ids:
                    allCOSeqs[i].append('-')
                else:
                    allCOSeqs[i].append(allDataSeqs[i][pos[i]])
                    pos[i] += 1
        else: # No insertion
            for i in range(numSeq):
                allCOSeqs[i].append(allDataSeqs[i][pos[i]])
                pos[i] += 1
        
        if all([pos[i] == len(allModelSeqs[i]) for i in range(numSeq)]):
            break
    
    if WRITEFILE:
        fname = os.path.splitext(os.path.basename(fname))[0]
        f = open(fname+".co", "w")
    for seq in allCOSeqs:
        print(''.join(seq))
        if WRITEFILE:
            f.write(''.join(seq)+'\n')
    if WRITEFILE:
        f.close()





def main():
    args = sys.argv[1:]
    fname = args[0]
    start_time = time.time()
    solve(fname)
    local_alignment_cpp(fname)
    end_time = time.time()
    print("#########################################################")
    print(f"Time Spent: {end_time-start_time} seconds")


if __name__ == "__main__":
    main()