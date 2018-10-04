from collections import defaultdict
from Bio import SeqIO
import numpy as np

class PrimerPair():
    def __init__(self,primerF,primerR,expectedPCR_size):
        self.primerF = primerF
        self.primerR = primerR
        self.expectedPCR_size = expectedPCR_size

    def getLineToFile(self):
        """
        :return: a line for output file
        """

        return "\t".join(self.primerF.getBasicInfo() + self.primerR.getBasicInfo()[1:] + [str(self.expectedPCR_size)])

    def getId(self):
        return self.primerF.id[:-1]

class Primer():
    """
    takes 1) fmt 6 file table from blast of primers vs genome. Use NCBI and download the results
    2) fasta of primers

    return three files:
    1. *single_PCR_primers.fasta file of primers that generate a single PCR product
    2. *._expectedPCR file where all verified PCR pairs and their products are written
    3. *single_PCR_primers.info where the information ONLY about primers generating single PCR product is presented
    """
    def __init__(self, chromosome, start, end, seq, id):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.seq = seq
        self.id = id
        self.primer_orientation = self.primer_orientation()

    def primer_orientation(self):
        side = 'forward'
        if self.start > self.end:
            side = 'reverse'
        return side

    def isForward(self):
        if self.primer_orientation == 'forward':
            return True
        else:
            return False

    def getEnd3(self):
        return self.end

    def getBasicInfo(self):
        return [self.chromosome, self.id,
                         str(self.seq.seq), str(self.start), str(self.end)]

class PrimerChecker():
    def __init__(self, blast6out, primerFasta):
        self.blast6out = blast6out
        self.primerFasta_name = primerFasta
        self.primerFasta = SeqIO.index(primerFasta, "fasta")
        self.mismatches = 2
        self.productSize = 10000
        self.PCR_pairs = [] # pairs of valid pcr primers as PrimerPair objects
        self.singleIDs = [] # ids of single PCR primer pairs
        self.main()


    def checkPrimer(self, sp_line):
        primer_length = len(self.primerFasta[sp_line[0]])
        #print(primer_length)
        if primer_length - self.mismatches > int(sp_line[3]):  # Alignment length
            print("Alignment length too small!")
            return False

        if self.mismatches < int(sp_line[4]): # NUMBER OF MISMATCHES
            print("Too many mismatches!")
            return False

        if int(sp_line[7]) != primer_length: # ckeck that 3' end is aligned
            print("3' end is not aligned correctly")
            return False
        return True

    def parseBlast(self):
        coordinates = {} #{primer_id:chromosome:Primer()}
        all_primers = []
        with open(self.blast6out) as infile:
            for lines in infile:
                if lines and not lines.startswith("#") and lines[0] != '':
                    #print(lines)
                    sp = lines.rstrip().split("\t")
                    all_primers.append(sp[0])
                    if self.checkPrimer(sp):
                        if sp[0] not in coordinates:
                            coordinates[sp[0]] = defaultdict(list)
                        coordinates[sp[0]][sp[1]].append(Primer(sp[1], int(sp[-4]), int(sp[-3]), self.primerFasta[sp[0]], sp[0]))

        print("Number of primers in BLAST file: ", len(set(all_primers)))
        print("Number of primers after BLAST: ", len(coordinates))
        return coordinates

    def getSinglePCRtab(self):
        with open(self.primerFasta_name + "single_PCR_primers.info", "w") as outTabSnigle:
            outTabSnigle.write("\t".join(["Chromosome", "F_primer_id", "F_primer_sequence", "F_Start", "F_end",
                                          "R_primer_id", "R_primer_sequence", "R_Start", "R_end", "Expected_PCR_product_length"]) + "\n")
            for pairs in self.PCR_pairs:
                if pairs.getId() in self.singleIDs:
                    outTabSnigle.write(pairs.getLineToFile() + "\n")

    def calculateClosestDistance(self,coordinates):
        d_product_per_primer_pair = defaultdict(int)
        with open(self.primerFasta_name + "_expectedPCR", "w") as outFile ,\
        open(self.primerFasta_name + "single_PCR_primers.fasta", "w") as singleOut:
            #print(coordinates)
            for primers in coordinates:
                #print(primers)
                #print(primers)
                if primers.endswith("F"):
                    paired_primer_id = primers[:-1] + "R"

                    #print(paired_primer_id)
                    if paired_primer_id in coordinates:
                        d_product_per_primer_pair[primers[:-1]] = 0

                        for chromosomes in coordinates[primers]:
                                if chromosomes in coordinates[paired_primer_id]:
                                    paire_pr_chromosome_positions = coordinates[paired_primer_id][chromosomes]
                                    #print(paire_pr_chromosome_positions[0].isForward())
                                    for primers_of_chromosome in coordinates[primers][chromosomes]:
                                        if primers_of_chromosome.isForward():
                                            primer_Set = [i for i in paire_pr_chromosome_positions if not i.isForward()]
                                            look_list_primer = [i.getEnd3() for i in primer_Set]
                                        else:
                                            primer_Set = [i for i in paire_pr_chromosome_positions if i.isForward()]
                                            look_list_primer = [i.getEnd3() for i in primer_Set]
                                            #print(look_list_primer)

                                        if look_list_primer:
                                            primer_Set = sorted(primer_Set, key=lambda x:x.getEnd3())
                                            look_list_primer = sorted(look_list_primer)
                                            index = np.searchsorted(look_list_primer, primers_of_chromosome.getEnd3())
                                            #print(index)
                                            #print(primers_of_chromosome.isForward())

                                            if not (index == 0 and primers_of_chromosome.isForward() == False) and \
                                                    not (index == len(look_list_primer) and primers_of_chromosome.isForward()):
                                                if not primers_of_chromosome.isForward():
                                                    closest_coordinate = index - 1
                                                else:
                                                    closest_coordinate = index

                                                paired_primer = primer_Set[closest_coordinate]
                                                original_primer = primers_of_chromosome


                                                #print(closest_coordinate)
                                                pcr_product = abs(primers_of_chromosome.getEnd3() - look_list_primer[closest_coordinate])
                                                #print(pcr_product)
                                                if pcr_product < self.productSize:
                                                    d_product_per_primer_pair[primers[:-1]] += 1
                                                    self.PCR_pairs.append(PrimerPair(original_primer, paired_primer,pcr_product))
                                                    outFile.write("\t".join(
                                                        [chromosomes,
                                                         primers,str(primers_of_chromosome.start), str(primers_of_chromosome.end),
                                                         paired_primer_id, str(primer_Set[closest_coordinate].start), str(primer_Set[closest_coordinate].end),
                                                         str(pcr_product)]
                                                    ) + "\n")
            print(d_product_per_primer_pair)
            cnt_many = []
            cnt_single = []
            for i, primers in enumerate(d_product_per_primer_pair):
                if d_product_per_primer_pair[primers] == 1:
                    print(i, primers)
                    print(d_product_per_primer_pair[primers])
                    print(">{0}F\n{1}\n".format(primers, str(self.primerFasta[primers + "F"].seq)))
                    print(">{0}R\n{1}\n".format(primers, str(self.primerFasta[primers + "R"].seq)))
                    cnt_single.append(primers)
                    singleOut.write(">{0}F\n{1}\n".format(primers, str(self.primerFasta[primers + "F"].seq)))
                    singleOut.write(">{0}R\n{1}\n".format(primers, str(self.primerFasta[primers + "R"].seq)))
                else:
                    cnt_many.append(primers)

            print("Total number of primer pairs with SINGLE pcr product:", len(cnt_single))
            print("Total number of primer pairs with MANY pcr products:", len(cnt_many))
            print("Total number of primer pairs WITHOUT detectable pcr product:", int(len(self.primerFasta)/2) - len(cnt_many) - len(cnt_single))
            self.singleIDs = cnt_single


    def main(self):
        coordinates = self.parseBlast()
        print(coordinates)
        self.calculateClosestDistance(coordinates)
        self.getSinglePCRtab()

if __name__ == "__main__":
    PrimerChecker(r'D:\Downloads\VBYV6YK8014-Alignment.txt', r'D:\PycharmProjects\MossRepeatomArticle\Primers_30.08.18_1.fasta')