import random as rand
import math
import numpy as np
import pandas as pd
from collections import Counter
import gzip
import os
import subprocess

import matplotlib.pyplot as plt
import seaborn as sns


class Futra:
    def __init__(self, sid, chrom,  start, strand, altchr, altpos, altstrand, svtype,  width, span, perSampleCount, inslen, gene1, gene2, qualityFilterPass, gap):
        self.sid = sid
        self.chrom = chrom
        self.start = start
        self.strand = strand
        self.altchr =  altchr
        self.altpos = altpos
        self.altstrand = altstrand
        self.svtype = svtype
        self.width = width
        self.span = span
        self.perSampleCount = perSampleCount
        self.inslen =  inslen
        self.gene1 = gene1
        self.gene2 = gene2
        self.qualityFilterPass = qualityFilterPass
        self.gap = gap
    def countSpace (self):
        if self.chrom == self.altchr:
            self.gap = abs(int(self.start) - int(self.altpos))
    def printVals(self):
        print ( "self.sid  = " + str(self.sid))
        print("self.chrom  = " + str(self.chrom))
        print("self.start  = " + str(self.start))
        print("self.strand  = " + str(self.strand))
        print("self.altchr  = " + str(self.altchr))
        print("self.altpos  = " + str(self.altpos))
        print("self.altstrand  = " + str(self.altstrand))
        print("self.svtype  = " + str(self.svtype))
        print("self.width  = " + str(self.width))
        print("self.span  = " + str(self.span))
        print("self.perSampleCount  = " + str(self.perSampleCount))
        print("self.inslen  = " + str(self.inslen))
        print("self.gene1  = " + str(self.gene1))
        print("self.gene2  = " + str(self.gene2))
        print("self.qualityFilterPass  = " + str(self.qualityFilterPass))



class BinThing:
    def __init__(self,idNumber_chrom, idNumber_total, chrom, start, end):
        self.idNumber_chrom = idNumber_chrom
        self.idNumber_total = idNumber_total
        self.idNumber_start = start
        self.end = end
        self.chrom = chrom


class TwoDimBinThing:
    def __init__(self, chrom1, chrom1Start, chrom1End, chrom2, chrom2Start, chrom2End, count, readsPerSampleTotal, futraCount, futras):
        self.futras = futras
        self.futraCount = futraCount
        self.chrom1 = chrom1
        self.chrom1Start =chrom1Start
        self.chrom1End =chrom1End
        self.chrom2 =chrom2
        self.chrom2Start = chrom2Start
        self.chrom2End = chrom2End
        self.count = count
        self.readsPerSampleTotal = readsPerSampleTotal






if __name__ == "__main__":
    print("start")
    filename = "C:/Users/Jacklyn/PycharmProjects/cse283Project/filtered_merged_10_30_2017.txt"

    binSize = 1000000
    lowerSepThreshold = 200000  # 200 kb go back and check tonight

    chrom9len = 138394717
    chrom9BinCount = int(chrom9len / binSize ) + 1

    bins = []
    currentPlace = 1
    for i in range(chrom9BinCount):
        print(i)
        bins.append([])
        start = currentPlace
        end = currentPlace + binSize
        currentPlace2 = 1
        for j in range(chrom9BinCount):
            start2 = currentPlace2
            end2 = currentPlace2 + binSize
            bin = TwoDimBinThing(9, start, end, 9, start2, end2, 0, 0, 0, [])
            bins[i].append(bin)
            currentPlace2 = currentPlace2 + binSize
        currentPlace = currentPlace + binSize



    print("chrom9 bins = " + str(chrom9BinCount))
    # chromLen = dict()
    # chromLen['1'] = 248956422
    # chromLen['2'] = 242193529
    # chromLen['3'] = 198295559
    # chromLen['4'] = 190214555
    # chromLen['5'] = 181538259
    # chromLen['6'] = 170805979
    # chromLen['7'] = 159345973
    # chromLen['8'] = 145138636
    # chromLen['9'] = 138394717
    # chromLen['10'] = 133797422
    # chromLen['11'] = 135086622
    # chromLen['12'] = 133275309
    # chromLen['13'] = 114364328
    # chromLen['14'] = 107043718
    # chromLen['15'] = 101991189
    # chromLen['16'] = 90338345
    # chromLen['17'] = 83257441
    # chromLen['18'] = 80373285
    # chromLen['19'] = 58617616
    # chromLen['20'] = 64444167
    # chromLen['21'] = 46709983
    # chromLen['22'] = 50818468
    #
    # chromLen['X'] = 156040895
    # chromLen['Y'] = 57227415
    # chromLen['M'] = 16569
    #
    # chromBinCount = dict()
    # for i in chromLen.keys():
    #     print(i)
    #     chromBinCount[i] = chromLen[i] % binSize
    #     rem = chromLen[i] - (chromBinCount[i] * binSize)
    #     if rem > 0:
    #         chromBinCount[i] = chromBinCount[i] + 1
    #
    # totalBinCount = 0
    #
    # for i in chromBinCount.keys():
    #     print(i)
    #     totalBinCount = totalBinCount + chromBinCount[i]
    #     print(totalBinCount)
    #
    # bins = []
    # totId = 1
    # iter = 0
    # for i in chromBinCount.keys():
    #     bincount = 1
    #     limit = chromLen[i]
    #     currentPlace = 1
    #     for j in range(chromBinCount[i]):
    #         chrom = i
    #         start = currentPlace
    #         end = currentPlace + binSize
    #         idNumber_chrom = bincount
    #         idNumber_total = totId
    #         bincount += 1
    #         totId += 1
    #         binNew = BinThing(idNumber_chrom, idNumber_total, chrom, start, end)
    #         currentPlace = end
    #         bins.append(binNew)
    #         iter += 1
    #
    # twoDimBins = [totalBinCount][totalBinCount]
    # for i in range(len(bins)):
    #     for j in range(len(bins)):
    #         twoDimBins[i][j] = TwoDimBinThing(bins[i].chrom, bins[i].start, bins[i].end, bins[j].chrom, bins[j].start,
    #                                           bins[j].end, 0, 0, 0, [])
    #
    with open(filename) as f:
        lines = f.readlines()
        print(lines)
    f.close()


    futras = []
    for i in lines:
        if len(i) > 0 and not i.startswith("sid"):
            cols = i.split("\t")
            print(cols)

            #    0                  1               2      3        4   5    6         7     8        9    10   11          12    13      14                        15
            #TCGA-05-4382	1030519:1_TCGA-05-4382	6	14490862	+	6	62330806	+	h2hINV	1030519	1	47839944	322		0		0						KHDRBS2	355394 bp to right of CD83 and 755908 bp to left of JARID2	5255627 bp to right of RAB23 and 60061 bp to left of KHDRBS2		FALSE	24	PASS	DSCRD	60	60	11	26.35	41	0	0	17	6:14490862(+)-6:62330806(+)__6_14479501_14504501D	LUAD	TRUE	TCGA	Primary	LUAD-Primary	LONG

            #  0                     1          2         3       4      5       6        7        8       9     10    11        12                  14      15     16          17     18       19       20              21             22            23       24      25      26    27     28            29      30      31          32      33      34      35      36       37     38
            #sid	                uidseqnames	chrom   start	strand	altchr	altpos	altstrand	svtype	ID	width	SPAN	per_sample_count	HOMSEQ	HOMLEN	INSERTION	INSLEN	gene1	gene2	cgc.gene1_100kb	cgc.gene2_100kb	gene1_100kb	gene2_100kb	break1	break2	fusion	L1	QUALITY	FILTER	EVDNC	DISC_MAPQ	MAPQ	TUMALT	TUMLOD	TUMDEP	NORMALT	NORMLOD	NORMDEP	SCTG	subtype	WGS	cohort	recurrent	subtype_recur	SPANclass

            #                 sid       chrom     start   strand  altchrom altstart  altstrand   svtype    width   SPAN   per sample count inslen    gene 1    gene 2  qualityFilterPass

            # newFutra = Futra("", 0, 0, "", 0, 0, "", "", "", "", 0, 0, "", "", False )
            # newFutra.sid = cols[0]
            # newFutra.chrom =cols[2]
            # newFutra.start =cols[3]
            # newFutra.strand =cols[4]
            # newFutra.altchr =cols[5]
            # newFutra.altpos =cols[6]
            # newFutra.altstrand =cols[7]
            # newFutra.svtype =cols[8]
            # newFutra.width =cols[10]
            # newFutra.span =
            # newFutra.perSampleCount =
            # newFutra.inslen =
            # newFutra.gene1 =
            # newFutra.gene2 =
            # newFutra.qualityFilterPass =


            newFutra = Futra(cols[0], cols[2], int(cols[3]), cols[4], cols[5], int(cols[6]), cols[7] ,  cols[8], cols[10], cols[11], cols[12], cols[17], cols[18], cols[19], cols[28], abs(int(cols[6]) - int(cols[3])))
            newFutra.countSpace()
            newFutra.printVals()
            futras.append(newFutra)






    for i in futras:
        if i.qualityFilterPass and i.gap <= lowerSepThreshold:

            if i.chrom is "9" and i.altchr is "9":
                for j in range(len(bins)):
                    for k in range(len(bins)):
                        binTemp = bins[j][k]

                        f1Match_b1 = False
                        f2Match_b1 = False
                        f2Match_b2 = False
                        f1Match_b2 = False
                        if i.start >= binTemp.chrom1Start and i.start <= binTemp.chrom1End:
                            f1Match_b1 = True
                        if i.altpos >= binTemp.chrom1Start and i.start <= binTemp.chrom1End:
                            f1Match_b2 = True
                        if i.start >= binTemp.chrom2Start and i.start <= binTemp.chrom2End:
                            f2Match_b1 = True
                        if i.altpos >= binTemp.chrom2Start and i.start <= binTemp.chrom2End:
                            f2Match_b2 = True

                        if f1Match_b1 and f2Match_b2:
                            binTemp.count += 1
                            binTemp.futras.append(i)
                            binTemp.readsPerSampleTotal = binTemp.readsPerSampleTotal + int(i.perSampleCount)
                        bins[j][k] = binTemp


    for i in range(len(bins)):
        for j in range(len(bins)):
            print(bins[i][j].count, end="\t")
        print("\n")

        contact_matrix = np.zeros((chrom9BinCount + 1, chrom9BinCount + 1))
        for i in range(len(bins)):
            for j in range(len(bins)):
                contact_matrix[i, j] = bins[i][j].count
        cell_type = "Tumor Fusions"

        ## This is Joaquin's code
    fig, ax = plt.subplots(figsize=(8, 6))

    vbounds = {'HFF': (0, 80), 'HEK': (0, 200), 'Fusion Pairs': (0, 50)}

    ax.set_title('Chromosome 9 Fusion Map', fontsize=14)
    sns.heatmap(contact_matrix + contact_matrix.T, cmap='bwr', vmin=vbounds[cell_type][0],vmax=vbounds[cell_type][1])


    fn = "C:/Users/Jacklyn/PycharmProjects/cse283Project_heat_map_fusions.png"
    fig.tight_layout()
    fig.savefig(fn, dpi=200)

