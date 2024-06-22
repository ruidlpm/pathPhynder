#!/usr/bin/env python3

# 14/05/2024

import argparse
import sys
import time
import re
import statistics
from collections import Counter


parser = argparse.ArgumentParser(description="estimates contamination from Ychr pileup")


#add options to argparser
parser.add_argument('-i', action="store", dest="pileup_input", type=str)
parser.add_argument('-o', action="store", dest="allele_count_output", type=str)
parser.add_argument('-d', action="store", dest="min_depth_cov", type=int)
parser.add_argument('-g', action="store", dest="isogglist", type=str)
parser.add_argument('-m', action="store", dest="mode", type=str)

#test parameters
try:
    options=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

pileup_input=options.pileup_input
allele_count_output=options.allele_count_output
min_depth_cov=options.min_depth_cov
mode=options.mode

if mode==None:
    mode='nofilter'
    print('mode: ' + mode)
elif mode!='nofilter' and mode!='transversions' :
    exit("-m parameter needs to be empty or -m 'nofilter' for considering all variants, or -m 'transversions' for considering non C/T G/A SNPs only.")
else:
    print('mode: ' + mode)

if min_depth_cov==None:
    min_depth_cov='1'
    print('min_depth_cov: ' + str(1))
else:
    print('min_depth_cov: ' + str(min_depth_cov))



def chunk_string(s, n):
    return [s[i:i+n] for i in range(len(s)-n+1)]


def get_bases(pileup):
    #first line is header - add hg(if any) and branch(later)
    parsed_pileup=[]
    for variant in pileup:
        CHR = variant[0]
        POS = variant[1]

        #clear '$' and '^*' patterns
        
        if int(variant[3]) ==0:
           variant.append('0')       
        #TODO coverage
        else:          
            while True:
                match = re.search(r"[+-](\d+)", variant[4])
                if match is None:
                    break
                variant[4] = variant[4][:match.start()] + variant[4][match.end() + int(match.group(1)):]
              
        variant[4] = variant[4].replace('$','')
        if '^' in variant[4]:
            variant[4] = variant[4].replace(variant[4][variant[4].index('^'):variant[4].index('^')+2],'')


        counter = Counter(chunk_string(variant[4], 1))
        REF_allele_count = counter[','] + counter['.']

        REF_allele = variant[2]

        REF_allele = REF_allele.upper()


        if REF_allele == 'A':
            A_count=REF_allele_count
        else:
            A_count=counter['A'] + counter['a']
        if REF_allele == 'T':
            T_count=REF_allele_count
        else:
            T_count=counter['T'] + counter['t']
        if REF_allele == 'C':
            C_count=REF_allele_count
        else:
            C_count=counter['C'] + counter['c']
        if REF_allele == 'G':
            G_count=REF_allele_count
        else:
            G_count=counter['G'] + counter['g']
        parsed_pileup.append([CHR, POS, REF_allele,A_count,T_count,C_count,G_count])
    return(parsed_pileup)




#read in haplogroup determination results
pileup=[]
with open(pileup_input, 'r') as in_pileup:
    for line in in_pileup.readlines():
        line=line.strip()
        col=line.split('\t')
        pileup.append(col)
    in_pileup.close()

base_calls = get_bases(pileup)


def contamination_at_base(minor_count, major_count):
    cont_est=minor_count*100/(minor_count + major_count)
    return(cont_est)

# call_haplogroups routine

if options.isogglist == 'none':
    print('please provide valid ISOGG SNP list, see data/ folder')
else:
    snps_hg_file=options.isogglist

    hg_pileup=[]
    with open(pileup_input, 'r') as in_hgs_pileup:
        for line in in_hgs_pileup.readlines():
            line=line.strip()
            col=line.split('\t')
            hg_pileup.append(col)
    in_hgs_pileup.close()

    hg_dict = {}
    with open(snps_hg_file) as hg_file:
        excluded_lines=0
        for line in hg_file:
            cols=line.split()
            pos=cols[2]
            if cols[3] in ['A','T','C','G'] and cols[4] in ['A','T','C','G']:
                hg_dict[pos] = cols
            else:
                excluded_lines=excluded_lines+1
        print('\n\n')
        print('\tinvalid markers excluded from input haplogroup list: ' + ' ' +  str(excluded_lines))
        print('\tvalid markers for haplogroup calling: ' + ' ' +  str(len(hg_dict)))
        print('\n\n')

    # line of the form ['CHR', 'POS', 'REF','A','T','C','G'], need to keep dict
    base_dict={}
    base_dict['A']=3
    base_dict['T']=4
    base_dict['C']=5
    base_dict['G']=6


    anc_snps=0
    der_snps=0
    NA_snps=0
    included_SNPs=0
    contamination_vect=[]
    with open(allele_count_output.replace('.intree.txt','.hgs.txt'),'w') as outfile:
        hdr = '\t'.join(['Mean_contamination','nSNPs_included','nSNPs_excluded','mode','min_depth_cov'])
        outfile.writelines(hdr)
        outfile.writelines('\n')
        for line in get_bases(hg_pileup):
            hg_pos=line[1]
            if hg_pos in hg_dict:
                hg_match=(hg_dict[hg_pos])
                pos=hg_match[2]
                marker=hg_match[0]
                hg=hg_match[1]
                ANC=hg_match[3]
                DER=hg_match[4]
                ANC_count=line[base_dict[ANC]]
                DER_count=line[base_dict[DER]]
                cov=(line[3]+line[4]+line[5]+line[6] )

                if cov>=min_depth_cov:
                    if mode=='transversions':
                        if ANC=='C' and DER=='T':
                            NA_snps=NA_snps+1
                            next
                        elif ANC=='T' and DER=='C':
                            NA_snps=NA_snps+1
                            next
                        elif ANC=='A' and DER=='G':
                            NA_snps=NA_snps+1
                            next
                        elif ANC=='G' and DER=='A':
                            NA_snps=NA_snps+1
                            next
                        else:
                            
                            if ANC_count==0 and DER_count==0:
                                NA_snps=NA_snps+1
                            elif ANC_count>0 and DER_count>0:
                                # determine minor
                                if ANC_count > DER_count:
                                    minor_count=DER_count
                                    major_count=ANC_count
                                    est=contamination_at_base(minor_count, major_count)
                                    contamination_vect.append(est)
                                    included_SNPs=included_SNPs+1
                                elif DER_count > ANC_count:
                                    minor_count=ANC_count
                                    major_count=DER_count
                                    est=contamination_at_base(minor_count, major_count)
                                    contamination_vect.append(est)
                                    included_SNPs=included_SNPs+1
                                elif DER_count == ANC_count:
                                    minor_count=ANC_count
                                    major_count=DER_count
                                    est=contamination_at_base(minor_count, major_count)
                                    contamination_vect.append(est)
                                    included_SNPs=included_SNPs+1
                            else:
                                if (ANC_count>0):
                                    contamination_vect.append(0)
                                    included_SNPs=included_SNPs+1
                                elif (DER_count>0):
                                    contamination_vect.append(0)
                                    included_SNPs=included_SNPs+1
                    else:
                        
                        if ANC_count==0 and DER_count==0:
                            NA_snps=NA_snps+1
                        elif ANC_count>0 and DER_count>0:
                            # determine minor
                            if ANC_count > DER_count:
                                minor_count=DER_count
                                major_count=ANC_count
                                est=contamination_at_base(minor_count, major_count)
                                contamination_vect.append(est)
                                included_SNPs=included_SNPs+1
                            elif DER_count > ANC_count:
                                minor_count=ANC_count
                                major_count=DER_count
                                est=contamination_at_base(minor_count, major_count)
                                contamination_vect.append(est)
                                included_SNPs=included_SNPs+1
                            elif DER_count == ANC_count:
                                minor_count=ANC_count
                                major_count=DER_count
                                est=contamination_at_base(minor_count, major_count)
                                contamination_vect.append(est)
                                included_SNPs=included_SNPs+1
                        else:
                            if (ANC_count>0):
                                contamination_vect.append(0)
                                included_SNPs=included_SNPs+1
                            elif (DER_count>0):
                                contamination_vect.append(0)
                                included_SNPs=included_SNPs+1
                else:
                    NA_snps=NA_snps+1                
        # print(str(round(statistics.mean(contamination_vect),4)))
        # print(str(included_SNPs))
        # print(str(NA_snps))        
        # print(mode)

        new_line=[str(round(statistics.mean(contamination_vect),4))+'%',str(included_SNPs),str(NA_snps), str(mode), str(min_depth_cov)]
        outfile.writelines('\t'.join(new_line))
        outfile.writelines('\n')
    outfile.close()

        # print('\nancestral state: ' + str(anc_snps))
        # print('\tderived state: ' + str(der_snps))
        # print('\texcluded: ' + str(NA_snps))
        # print("\twritten to " + allele_count_output  + '\n\n')

        # print(round(statistics.stdev(contamination_vect),4 ))
