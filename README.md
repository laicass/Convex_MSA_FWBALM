# MSA Problem FW-BALM Solver
Solver for MSA problem using FW-BALM algorithm.
## How to run
```py MSA_Convex.py test\74.msa```
## Util functions
- Calculate Sum-of-Pair Score
```py util/co2SPScore.py [MSA file].co```
- Calculate Star Alignment Score
```py util/co2SPScore.py [MSA file].co```
- Plot FW steps comparison
```py util/plot_log.py [FW-ALM_log].log [FW-BALM_log].log```
## Output files
- .msa: MSA problem input
- .co : Clustal Omega formatted output
- .rec: Recovered string
- .log: Log output
