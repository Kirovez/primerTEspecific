# primerTEspecific
python scripts to design specific primers for transposons


The program was developed in frame of the RNF grant (№17-14-01189)

Our website: [Plantprotlab.com](http://plantprotlab.com/)
## Getting Started

### Prerequisites

primerTEspecific requires:
* python v3.6
python packages to be installed: biopython, numpy
(run command: `pip install biopython networkx`)

### Run primerTEspecific

`python checkPrimerSpecificity.py [-h] fmt6blast primers`

* fmt6blast - results of blastn (format 6 table) primers vs genome
* primers - fasta file with primers. Note! Two primers must differ in single end
              letter: "F" for one and "R" for another primer