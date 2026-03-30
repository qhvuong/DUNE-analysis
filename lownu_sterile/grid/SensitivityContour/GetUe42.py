import sys
import numpy as np


if __name__ == "__main__":

    n = int(sys.argv[1])
    N = int(sys.argv[2])
    Ue42_min = float(sys.argv[3])
    Ue42_max = float(sys.argv[4])
    Um42_min = float(sys.argv[5])
    Um42_max = float(sys.argv[6])
    Ut42     = float(sys.argv[7])
    dm2_min  = float(sys.argv[8])
    dm2_max  = float(sys.argv[9])

    log_min = np.log10(Ue42_min)
    log_max = np.log10(Ue42_max)
    binWidth = (log_max - log_min) / N
    
    i = int(n%N)
    
    #bin_edge_lower = 10 ** (log_min + i * binWidth)
    #bin_edge_upper = 10 ** (log_min + (i+1) * binWidth)
    
    bin_center = 10 ** (log_min + (i+0.5) * binWidth)

    print(bin_center)

