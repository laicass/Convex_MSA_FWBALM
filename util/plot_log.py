import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    log_files = sys.argv[1:]
    fig, ax1 = plt.subplots()
    #ax2 = ax1.twinx()
    for f, label in zip(log_files, ["FW-ALM", "FW-BALM"]):
        log = np.loadtxt(f)
        CoZ_log = log[:,-1]
        fw_iter = np.sum(log[:,:-1], axis=1)
        total_step = int(np.sum(fw_iter))
        ax1.bar(np.arange(len(fw_iter)), fw_iter, label=label+" total steps:"+str(total_step))
        #ax2.plot(CoZ_log)

    plt.title("Number of FW step comparison")
    plt.xlabel("Outer loop iteration")
    ax1.set_ylabel("FW steps per iteration")
    #ax2.set_ylabel("Objective value", color='g')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()