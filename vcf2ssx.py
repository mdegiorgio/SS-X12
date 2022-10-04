#Converts a .vcf.gz to .txt.gz, which will be read into SS-X12
#If your file is in VCF format, USE THIS SCRIPT! All other formats, convert to VCF first!
#Check the SS-X12 manual for operation instructions
#You will need to provide the indices of important summary data within the VCF data lines...if you're not formatted exactly as VCFv4.2
#Outfile name will be the same as infile name, but with different file extension (.vcf or .vcf.gz is converted to .txt.gz)

#This script is a bit slow...the 6M+ lines of 1000G chr1 took 9.5 hours to complete, so allow appropriate time for vcf2ssx to run

#Consider removing sites from your VCF that are NOT biallelic SNPs...vcf2ssx will run MUCH faster

import sys
import subprocess
import gzip
import re

def conversion(inline, chromID, pos, snpID, ref, alt, fdat): #current data line, identity of chromosome (index), physical position of SNP (index), rsID of SNP (index), reference allele (index), alternate allele (index), first data position (index)
    reformat = []
    location = 0
    for item in inline:
        if ((location == chromID) or (location == pos) or (location == ref) or (location == alt)):
            reformat.append(item)
        elif location == snpID:
            reformat.insert(1, item)
        elif location >= fdat:
            if (("0|0" in item) or ("0/0" in item)): #homozygous reference
                reformat.append("1")
            elif "0|1" in item: #heterozygous "M"
                reformat.append("2")
            elif "1|0" in item: #heterozygous "P"
                reformat.append("3")
            elif (("1|1" in item) or ("1/1" in item)): #homozygous alternate
                reformat.append("4")
            elif (("0/1" in item) or ("1/0" in item)): #heterozygous unphased
                reformat.append("5")
            elif "." in item: #missing data encoded as a dot in VCF 4.2 (and presumably others)
                reformat.append("N")

        location += 1

    return(reformat)

invcf_name = sys.argv[1]

if ".gz" in invcf_name:
    invcf = gzip.open(invcf_name, "rb") #.vcf.gz file that we are reading
else:
    invcf = open(invcf_name, "r+")

orthodox = sys.argv[2] #either yes or no...if it's orthodox, we're not requiring more arguments; yes if the header is "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT"

if orthodox == "no": #give the indices (starting from 0 for the first) in the data lines that each part of the summary is found
    chromosome = int(sys.argv[4])
    position = int(sys.argv[5])
    rsID = int(sys.argv[6])
    reference = int(sys.argv[7])
    alternate = int(sys.argv[8])
    first_data = int(sys.argv[9])
else: #standard arrangement of indices
    chromosome = 0
    position = 1
    rsID = 2
    reference = 3
    alternate = 4
    first_data = 9

outname = invcf_name[:invcf_name.index(".vcf")] + ".txt"
outdest = open(outname, "a+")
sys.stdout = outdest

thereyet = 0
while thereyet != "#CHROM": #get to the header line
    thereyet = (invcf.readline().split())[0]

nonbi_removed = sys.argv[3]

if nonbi_removed == "no":
    while True:
        theline = invcf.readline().split()

        vartype = map(lambda x: re.search(r"VT=SNP", x), theline) #search for a match in theline

        try:
            if type(theline[0]) == str:
                if len(set(vartype)) > 1: #if it's 0, then there wasn't a match
                    macheck = map(lambda x: re.search(r"MULTI_ALLELIC", x), theline)
                    if len(set(macheck)) == 1:
                        converted = conversion(theline, chromosome, position, rsID, reference, alternate, first_data)
                        print(" ".join(entry for entry in converted))
        except:
            sys.stdout = sys.__stdout__
            message = "Analysis of data file \'" + invcf_name + "\' is complete!"
            print(message)
            break
else: #no need to do regular expression matching if all of the input lines are of the correct type
    while True:
        theline = invcf.readline().split()

        try:
            if type(theline[0]) == str:
                converted = conversion(theline, chromosome, position, rsID, reference, alternate, first_data)
                print(" ".join(entry for entry in converted))
        except:
            sys.stdout = sys.__stdout__
            message = "Analysis of data file \'" + invcf_name + "\' is complete!"
            print(message)
            break

outdest.close()

invcf.close()

zipcom = "bgzip " + outname

subprocess.call(zipcom, shell=True)
