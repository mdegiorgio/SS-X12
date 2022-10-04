# Welcome to SS-X12! This is a script written in Python 2.7 that searches for signatures of shared selective sweeps in two or more populations, from haplotype, OR multilocus genotype data (though it works better for haplotypes)
# Please consult the user manual included with this software for operation instructions including command line and input file formatting, as well as example data
# I'm sure there is plenty that can still be optimized/improved, so please let me know if you'd like to see a particular update or type of update (or if you program in a faster language than Python and can help make this script run more quickly)
# The phrase "sequence strings" used throughout the comments refers to the haplotype *OR* MLG (it's a generic term to reflect that either input is fine)
# References and variable assignments to H12 and H123 are equivalent to G12 and G123 if input data are MLG-formatted

# This script was written by Alexandre M. Harris (contact: amh522@psu.edu OR xistheway@gmail.com), please cite appropriately if you use this script
# SS-X12 is free to use and modify (once again, with the proper attribution)

import sys
import numpy as np
from collections import defaultdict
import gzip

def H1calc(x): #x is a list of frequencies
    sqlist = []
    for item in x:
        item = item ** 2
        sqlist.append(item)
    return(sum(sqlist)) #expected homozygosity

def H2calc(x):
    sqlist = []
    for item in x:
        item = item ** 2
        sqlist.append(item)
    return(sum(sqlist) - max(sqlist)) #minus p1^2

def H12calc(x):
    sqlist = []
    x = sorted(x)
    for item in x:
        item = item ** 2
        sqlist.append(item)
    if len(x) == 1:
        return(1)
    else:
        return(sum(sqlist) + (2 * x[len(x) - 1] * x[len(x) - 2])) #plus 2p1p2

def H123calc(x):
    sqlist = []
    x = sorted(x)
    for item in x:
        item = item ** 2
        sqlist.append(item)
    if len(x) == 1:
        return(1)
    elif len(x) == 2:
        return(sum(sqlist) + (2 * x[len(x) - 1] * x[len(x) - 2])) #if there aren't enough frequencies, then reduces to H12
    else:
        return(sum(sqlist) + (2 * x[len(x) - 1] * x[len(x) - 2]) + (2 * x[len(x) - 3] * x[len(x) - 2]) + (2 * x[len(x) - 1] * x[len(x) - 3]))  #plus 2p1p2, 2p1p3, 2p2p3

def dictfreqs(instrings, breakpoint, gpF, gpS): #sequence strings, population breaks, first group, second group; obtains frequencies which will be used to compute our statistics
    countdict = defaultdict(int) #total count in population
    countdict_popk = {} #dict of countdicts for each population [each entry formatted like above countdict]

    for nps in range(len(breakpoint)): #an entry into countdict for each population
        countdict_popk["countdict_pop" + str(nps + 1)] = defaultdict(int) #each entry in the dict is indexed as 'countdict_pop1/2/3/4...'

    c_index = 0 #helps me decide where an entry goes into countdict_popk

    for string in instrings:
        countdict[string] += 1 #increase count of that sequence string

        whichpop = breakpoint.index(int(min([bp for bp in breakpoint if bp > c_index]))) + 1 #which population is this, based on the current line?

        for cd in countdict_popk: #each cd is the key (a string) corresponding to a value in the dict
            if cd == "countdict_pop" + str(whichpop):
                countdict_popk[cd][string] += 1 #increment the count in that population
            else:
                countdict_popk[cd][string] += 0 #acknowledge the existence of the sequence string (haplo/mlg) in the pooled population

        c_index += 1

    datalen = float(len(instrings)) #total pooled sample size

    datafreqs = [] #frequencies of sequence string within the window for pooled sample (to calc [H1, H2, H12, H123]tot)

    datafreqs_SCT = {} #each entry is the datafreqs of each pair of populations (so each entry is a list like above, used to get the H12tot for the pair); note that the order of sequence strings is not necessarily the same across dict elements (but within an element, we're definitely averaging the frequency of the same sequence string)

    datafreqs_sqdiff_SCT = {} #each entry is a list of sqdiffs for each sequence string for each pair of populations (used to get the fdiff for the pair)

    datafreqs_eachpop = {} #we will ultimately want the H12 for each population to get our min/max ratio for SS-H12

    if ((gpF != "none") and (gpS != "none")): #only exists if we're doing H_diff based on grouping
        datafreqs_sqdiff_GROUPS = []
        datafreqs_eachgroup = {} #we want to get H12(3) for both groups from the datafreqs so that we can correct H12anc

        datafreqs_eachgroup["groupF"] = []
        datafreqs_eachgroup["groupS"] = []

        gpF = gpF.split(",")
        gpS = gpS.split(",")

        gpFpops = [] #names of the keys in countdict_popk belonging to group 1
        gpSpops = []

        for gf in gpF: #fill the frequencies for first group
            cdname = "countdict_pop" + gf
            gpFpops.append(cdname)

        for gs in gpS:
            cdname = "countdict_pop" + gs
            gpSpops.append(cdname)

    for item in countdict.values(): #get frequencies for pooled population
        freq = float(item) / datalen
        datafreqs.append(freq)

    for itepo in countdict_popk["countdict_pop1"]: #for each sequence string in the first population [note that all sequence string counts are included in each list, even if zero, so we're interrogating all of them]
        freqouts = {} #all the frequencies of that sequence string across each subpopulation; works because we're defining a count for each haplo regardless of its presence (zero if not there)
        freqsgroupF = [] #frequency of itepo in first group [list size = number of pops in first group]
        freqsgroupS = []

        sampsizesF = [] #sample sizes of each constituent population of the first group (for weighting)
        sampsizesS = []

        for cd in countdict_popk: #for each population in the nested dict
            frepk = float(countdict_popk[cd][itepo]) / sum(countdict_popk[cd].values()) #sum of the values here is total count of sequence strings for that population
            freqouts[cd] = frepk #the frequency of sequence string itepo [the current one we're interrogating] in population cd is added to freqouts
            if ((gpF != "none") and (gpS != "none")):
                if cd in gpFpops: #we can confirm whether the character string that is a key in the nested dict appears as an element of one of the gp lists
                    freqsgroupF.append(frepk) #adds that frequency to the first group if the current pop is in it
                    sampsizesF.append(float(sum(countdict_popk[cd].values())))
                elif cd in gpSpops:
                    freqsgroupS.append(frepk)
                    sampsizesS.append(float(sum(countdict_popk[cd].values())))

        for itemF in freqouts: #keys (populations) in freqouts
            if itemF not in datafreqs_eachpop.keys():
                datafreqs_eachpop[itemF] = [] #create empty list for sequence string frequencies in the individual pop
            datafreqs_eachpop[itemF].append(freqouts[itemF]) #add the frequency of this sequence string to the growing list
            for itemS in freqouts: #keys in freqouts (again--reiterating)
                if itemF != itemS: #they can't be identical, but we also can't have repeats
                    named = itemF + "_" + itemS #intended key for this value
                    altnamed = itemS + "_" + itemF #key that may already exist and therefore contain this value
                    if altnamed not in datafreqs_sqdiff_SCT.keys():
                        sqdiff = (freqouts[itemF] - freqouts[itemS]) ** 2 #squared difference in the frequency of this sequence string in popF and popS (fdiff)
                        if named not in datafreqs_sqdiff_SCT.keys(): #need to have an initial condition
                            datafreqs_sqdiff_SCT[named] = [] #create the list to which sqdiffs append
                        datafreqs_sqdiff_SCT[named].append(sqdiff)
                    if altnamed not in datafreqs_SCT.keys():
                        totalsize = sum(countdict_popk[itemF].values()) + sum(countdict_popk[itemS].values())
                        weightF = float(sum(countdict_popk[itemF].values())) / float(totalsize)
                        weightS = float(sum(countdict_popk[itemS].values())) / float(totalsize)
                        wavg = (freqouts[itemF] * weightF) + (freqouts[itemS] * weightS) #weighted average accounting for sample sizes
                        if named not in datafreqs_SCT.keys():
                            datafreqs_SCT[named] = []
                        datafreqs_SCT[named].append(wavg)

        if ((gpF != "none") and (gpS != "none")):
            poolsizeF = float(sum(sampsizesF))
            poolsizeS = float(sum(sampsizesS))

            wavgFG = []
            wavgSG = []

            for ggrF in range(len(freqsgroupF)):
                weightF = sampsizesF[ggrF] / poolsizeF
                weightfreqF = freqsgroupF[ggrF] * weightF
                wavgFG.append(weightfreqF)

            for ggrS in range(len(freqsgroupS)):
                weightS = sampsizesS[ggrS] / poolsizeS
                weightfreqS = freqsgroupS[ggrS] * weightS
                wavgSG.append(weightfreqS)

            wavgFG = sum(wavgFG)
            wavgSG = sum(wavgSG)

            sqdiffgroups = (wavgFG - wavgSG) ** 2 #squared difference in the frequency of the sequence string in the two lineages
            datafreqs_sqdiff_GROUPS.append(sqdiffgroups)

            datafreqs_eachgroup["groupF"].append(wavgFG)

            datafreqs_eachgroup["groupS"].append(wavgSG)

    if ((gpF != "none") and (gpS != "none")):
        return(datafreqs, datafreqs_SCT, datafreqs_sqdiff_SCT, datafreqs_eachpop, datafreqs_sqdiff_GROUPS, datafreqs_eachgroup)
    else:
        return(datafreqs, datafreqs_SCT, datafreqs_sqdiff_SCT, datafreqs_eachpop)

def get_freqs(matr, breakpoint, gpF, gpS, ploidy): #data matrix (individuals x SNPs), indices for start of next population, first group, second group; turns your input data into strings interpretable by the script [[CHANGE THIS SO THAT IT WORKS WITH HAPLO (1/2/3/4/N) OR MLG (1/4/5)]]
    strings = []
    for ind in range(matr.shape[0]): #iterate through rows, which after transpose are individuals
        if ploidy == "hap":
            current_string_h1 = 0 #variable for building first haplotype string from example haplotype data (must be reset for each individual)
            current_string_h2 = 0
            for item in matr[ind]: #SNPs associated with that individual
                if current_string_h1 == 0: #first loop (no SNPs added for that individual yet)
                    if item == "1":
                        current_string_h1 = "0"
                        current_string_h2 = "0"
                    elif item == "2":
                        current_string_h1 = "0"
                        current_string_h2 = "1"
                    elif item == "3":
                        current_string_h1 = "1"
                        current_string_h2 = "0"
                    elif item == "4":
                        current_string_h1 = "1"
                        current_string_h2 = "1"
                    elif item == "N": #missing data identifier--change this to whatever your missing data identifier is! [default behavior for missing data is to generate a new sequence string...this conservatively reduces the occurrence of outlying values for windows with missing sites]
                        current_string_h1 = item
                        current_string_h2 = item
                else: #concatenate current SNP
                    if item == "1":
                        current_string_h1 = current_string_h1 + "0"
                        current_string_h2 = current_string_h2 + "0"
                    elif item == "2":
                        current_string_h1 = current_string_h1 + "0"
                        current_string_h2 = current_string_h2 + "1"
                    elif item == "3":
                        current_string_h1 = current_string_h1 + "1"
                        current_string_h2 = current_string_h2 + "0"
                    elif item == "4":
                        current_string_h1 = current_string_h1 + "1"
                        current_string_h2 = current_string_h2 + "1"
                    elif item == "N":
                        current_string_h1 = current_string_h1 + item
                        current_string_h2 = current_string_h2 + item

                if len(current_string_h2) == len(matr[ind]): #compare the number of characters (SNPs) in string sequence to length of row in array to indicate completion of [haplotype]
                    strings.append(current_string_h1)
                    strings.append(current_string_h2)
        elif ploidy == "mlg":
            current_string_mlg = 0
            for item in matr[ind]:
                if current_string_mlg == 0:
                    current_string_mlg = item
                else:
                    current_string_mlg = current_string_mlg + item

                if len(current_string_mlg) == len(matr[ind]):
                    strings.append(current_string_mlg)

    return(dictfreqs(strings, breakpoint, gpF, gpS)) #tuple; output of dictfreqs is datafreqs list and sqdiff list

def update_SNPs(matchsites, splitpop, VALall, INFOall, SNPall, SNPcurr): #advances the analysis window
    popsnps = [] #un-concatenated list of allelic states in the sample at current SNP [example file: 1-4]

    for ind in matchsites: #gets valid indices for chosen populations
        popsnps.append(splitpop[ind + 5]) #data starts at index 5 in SS-X12 file

    popsnps_set = set(popsnps) #number of site types
    if ((len(popsnps_set) > 1) or (("2" or "3") in popsnps) or ("5" in popsnps)): #need to consider haplotype polymorphism, so keeping 2 and 3 [change this as necessary, or remove, to fit the format of your data...in the example, we don't want to include invariant sites, which would be those consisting of only individuals encoded as 1 or 4]
        VALall.append(popsnps) #for this SNP, report the allelic states
        INFOall.append(SNPcurr) #appends if informative in study population
    SNPall.append(SNPcurr) #appends regardless of informativeness

def calc_stats(popdata, outname, wincent, SNPall, INFOall, breakpoint, gpF, gpS, ploidy): #computes all statistics
    t_current_values = np.transpose(np.array(popdata)) #becomes an array of individuals x SNPs

    all_calc = get_freqs(t_current_values, breakpoint, gpF, gpS, ploidy) #spits out a tuple of lists

    whereto = outname + ".txt" #change as necessary

    try:
        fr_total = all_calc[0] #list of all freqs for pool
        fr_SCT = all_calc[1] #dict of datafreqs for each pair of populations
        fr_diff_SCT = all_calc[2] #dict of sqdiffs for each pair of populations
        fr_eachpop = all_calc[3] #dict of datafreqs for each pop in the analysis

        if ((gpF != "none") and (gpS != "none")):
            fr_diff_group = all_calc[4] #list of diffs between the two groups
            fr_eachgroup = all_calc[5] #dict of datafreqs for both groups

        Hdiff_SCT = {}
        #######################################
        H12tot_SCT = {}
        H123tot_SCT = {}
        H12anc_SCT = {}
        H123anc_SCT = {}
        #######################################
        H12_eachpop = {}
        H123_eachpop = {}
        H12_eachpop_ratio = {}
        H123_eachpop_ratio = {}
        H12_anc_corr = {}
        H123_anc_corr = {}

        if (fr_total != []): #the window contains at least an informative site
            H12total_all = H12calc(fr_total)
            H123total_all = H123calc(fr_total)
            H2H1total_all = H2calc(fr_total) / H1calc(fr_total)
        else:
            H12total_all = 1.0 #no diversity yields complete homozygosity
            H123total_all = 1.0
            H2H1total_all = 0.0

        if ((gpF != "none") and (gpS != "none")):
            H12_eachgroup = {}
            H123_eachgroup = {}

            H12_eachgroup["groupF"] = H12calc(fr_eachgroup["groupF"])
            H12_eachgroup["groupS"] = H12calc(fr_eachgroup["groupS"])
            H123_eachgroup["groupF"] = H123calc(fr_eachgroup["groupF"])
            H123_eachgroup["groupS"] = H123calc(fr_eachgroup["groupS"])

            Hdiff_group = sum(fr_diff_group)

            H12_anc_group = H12total_all - Hdiff_group
            H12_anc_group_ratio = min(H12_eachgroup.values()) / max(H12_eachgroup.values())
            H12_anc_group_corr = H12_anc_group * H12_anc_group_ratio #SS-H12
            #######################################
            H123_anc_group = H123total_all - Hdiff_group
            H123_anc_group_ratio = min(H123_eachgroup.values()) / max(H123_eachgroup.values())
            H123_anc_group_corr = H123_anc_group * H123_anc_group_ratio  #SS-H123

        for sqd in fr_diff_SCT:
            Hdiff_SCT[sqd] = sum(fr_diff_SCT[sqd])

        for dfs in fr_SCT:
            if len(fr_SCT[dfs]) > 1:
                H12total_popjk = H12calc(fr_SCT[dfs])
                H123total_popjk = H123calc(fr_SCT[dfs])
            else:
                H12total_popjk = 1.0
                H123total_popjk = 1.0

            H12tot_SCT[dfs] = H12total_popjk
            H123tot_SCT[dfs] = H123total_popjk
            H12anc_thispair = H12total_popjk - Hdiff_SCT[dfs]
            H123anc_thispair = H123total_popjk - Hdiff_SCT[dfs]
            H12anc_SCT[dfs] = H12anc_thispair
            H123anc_SCT[dfs] = H123anc_thispair

        for dfesF in fr_eachpop:
            for dfesS in fr_eachpop:
                if dfesF != dfesS:
                    namedd = dfesF + "_" + dfesS
                    altnamedd = dfesS + "_" + dfesF
                    if altnamedd not in H12_eachpop_ratio.keys():
                        H12_thatpopF = H12calc(fr_eachpop[dfesF])
                        H12_thatpopS = H12calc(fr_eachpop[dfesS])
                        #######################################
                        H123_thatpopF = H123calc(fr_eachpop[dfesF])
                        H123_thatpopS = H123calc(fr_eachpop[dfesS])

                        H12_eachpop[dfesF] = H12_thatpopF
                        H12_eachpop[dfesS] = H12_thatpopS
                        #######################################
                        H123_eachpop[dfesF] = H123_thatpopF
                        H123_eachpop[dfesS] = H123_thatpopS

                        H12_eachpop_ratio[namedd] = min(H12_thatpopF,H12_thatpopS) / max(H12_thatpopF,H12_thatpopS)
                        H12_anc_corr[namedd] = H12anc_SCT[namedd] * H12_eachpop_ratio[namedd]
                        #######################################
                        H123_eachpop_ratio[namedd] = min(H123_thatpopF,H123_thatpopS) / max(H123_thatpopF,H123_thatpopS)
                        H123_anc_corr[namedd] = H123anc_SCT[namedd] * H123_eachpop_ratio[namedd]

        c_H12anc_max = max(H12anc_SCT.values())
        c_H12anc_corr_max = max(H12_anc_corr.values())
        #######################################
        c_H123anc_max = max(H123anc_SCT.values())
        c_H123anc_corr_max = max(H123_anc_corr.values())

        c_H12anc_min = min(H12anc_SCT.values())
        c_H12anc_corr_min = min(H12_anc_corr.values())
        #######################################
        c_H123anc_min = min(H123anc_SCT.values())
        c_H123anc_corr_min = min(H123_anc_corr.values())

        if ((gpF != "none") and (gpS != "none")):
            output = [H12total_all, H123total_all, H2H1total_all, c_H12anc_corr_max, c_H123anc_corr_max, c_H12anc_corr_min, c_H123anc_corr_min, H12_anc_group_corr, H123_anc_group_corr, wincent, len(SNPall), len(INFOall)]
        else:
            output = [H12total_all, H123total_all, H2H1total_all, c_H12anc_corr_max, c_H123anc_corr_max, c_H12anc_corr_min, c_H123anc_corr_min, wincent, len(SNPall), len(INFOall)]
    except:
        if ((gpF != "none") and (gpS != "none")):
            output = [1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, wincent, len(SNPall), len(INFOall)]
        else:
            output = [1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, wincent, len(SNPall), len(INFOall)]

    sendto = open(whereto, "a+")
    sys.stdout = sendto
    print(" ".join(str(entry) for entry in output))
    sendto.close()

invcf_name = sys.argv[1]

if ".gz" in invcf_name:
    chromfile = gzip.open(invcf_name, "rb") #name the properly-formatted infile to read...name needs to be formatted as _hap.txt or _mlg.txt (with no other periods)
else:
    chromfile = open(invcf_name, "r+")

popselect = sys.argv[2] #name specific population(s) from the second line
headfile = open(sys.argv[3], "r+") #header file containing the population assignment for each individual in the data
outdest = sys.argv[4] #just include chromosome and date (no file extension given by user!)
size_window = int(sys.argv[5]) #how big is the window (kb)?
shift_window = int(sys.argv[6]) #shift window by how many kb?
groupF = sys.argv[7] #If we're grouping, who's in the first group? (If there are 3 groups [GIH, CEU, YRI], then the group IDs are 1, 2, 3 in the order that they're specified)
groupS = sys.argv[8]

eachpop = popselect.split(",")

populations = (headfile.readline()).split()

poparray = np.array(populations)

okinds = [] #nested list of indices of each individual in each population sample, formatted as [[pop1_ind1, pop1_ind2,...],...[popX_ind1,...]]

whichploidy = invcf_name[invcf_name.index(".") - 3:invcf_name.index(".")] #will either be "hap" or "mlg"

bpoints = [] #this is what we're using because the example file consists of haplotype data

for po in eachpop:
    match = np.flatnonzero(poparray == po)
    match = match.tolist() #indices for that population
    okinds.append(match)

    if len(bpoints) != 0: #if this isn't the first breakpoint
        if whichploidy == "hap":
            bpoints.append((len(match)) * 2 + bpoints[len(bpoints) - 1]) #haplotype breakpoints
        elif whichploidy == "mlg":
            bpoints.append(len(match) + bpoints[len(bpoints) - 1]) #individual breakpoints
    else: #if this is the first breakpoint
        if whichploidy == "hap":
            bpoints.append(len(match) * 2)
        elif whichploidy == "mlg":
            bpoints.append(len(match))

okinds = [item for sublist in okinds for item in sublist] #flatten okinds

firstloc_position = chromfile.tell() #this should just be zero, but it's a legacy command

first_SNP = int((chromfile.readline()).split()[2]) #saves the position of that SNP *in the chromosome* [needed for window-building]

chromfile.seek(firstloc_position) #returns pointer to 1st line so that the first SNP isn't omitted in the loop

window_center = first_SNP + (size_window / 2) #we are using a nucleotide-delimited sliding window approach, and the first window's lower bound is the first SNP

current_values = [] #list of lists (sublists are genotypes at that SNP) [co-indexed with current_infosites] [size S x n for S SNPs across n individuals]

current_SNPs = [] #SNPs in a window (reset every time; list of positions, whether or not they're informative in the study populations)

current_infosites = [] #informative sites in a window (reset every time; list of positions, but only polymorphic in study populations)

next_winstart = 0 #appends the index of the first valid element after shift
next_infowinstart = 0 #same, but for informative site [physical position of the window is the same for both this and next_winstart, but next_infowinstart is likely going to be a subset of next_winstart]

while True:
    position = chromfile.readline() #starts back at first SNP

    window_size = range((window_center - (size_window / 2)), (window_center + ((size_window / 2) + 1))) #physical definition of the window

    if (position != ''): #if a new line (SNP) has been read (we're not at the end of the file)
        posplit = position.split()
        current_SNP = int(posplit[2]) #this indexing works for the example file, but change it to match your line formatting

        if (current_SNP in window_size):
            update_SNPs(okinds, posplit, current_values, current_infosites, current_SNPs, current_SNP)

        else:
            calc_stats(current_values, outdest, window_center, current_SNPs, current_infosites, bpoints, groupF, groupS, whichploidy)

            update_SNPs(okinds, posplit, current_values, current_infosites, current_SNPs, current_SNP)

            for item in current_SNPs:
                if int(item) >= (min(window_size) + shift_window): #next valid winstart has to be at least one winshift from previous winstart
                    next_winstart = current_SNPs.index(item)
                    break #don't continue after finding the next valid SNP

            for iitem in current_infosites:
                if int(iitem) >= (min(window_size) + shift_window):
                    next_infowinstart = current_infosites.index(iitem)
                    break

            window_center += shift_window
            window_size = range((window_center - (size_window / 2)), (window_center + ((size_window / 2) + 1))) #needed to shift window

            current_SNPs = current_SNPs[next_winstart:]
            current_infosites = current_infosites[next_infowinstart:]
            current_values = current_values[next_infowinstart:]

            while (current_SNP not in window_size): #scans and reports until a new window contains the current SNP
                a_current_values = current_values[:len(current_values) - 1] #doesn't include non-fitting SNP

                calc_stats(a_current_values, outdest, window_center, current_SNPs, current_infosites, bpoints, groupF, groupS, whichploidy)

                for item in current_SNPs:
                    if int(item) >= (min(window_size) + shift_window):
                        next_winstart = current_SNPs.index(item)
                        break

                for iitem in current_infosites:
                    if int(iitem) >= (min(window_size) + shift_window):
                        next_infowinstart = current_infosites.index(iitem)
                        break

                window_center += shift_window
                window_size = range((window_center - (size_window / 2)), (window_center + ((size_window / 2) + 1)))

                current_SNPs = current_SNPs[next_winstart:]
                current_infosites = current_infosites[next_infowinstart:]
                current_values = current_values[next_infowinstart:]
    else: #likely to create an incomplete window...filter this out if it's too small
        calc_stats(current_values, outdest, window_center, current_SNPs, current_infosites, bpoints, groupF, groupS, whichploidy)

        sys.stdout = sys.__stdout__
        message = "Analysis of data file \'" + invcf_name + "\' is complete!"
        print(message)
        break

chromfile.close()
