#!/usr/bin/env python3

import argparse
import sys
import time
import re
from collections import Counter

parser = argparse.ArgumentParser(description="Parses samtools pileup, accounts for deamination and mismatches, outputs allele count")


#add options to argparser
parser.add_argument('-i', action="store", dest="pileup_input", type=str)
parser.add_argument('-m', action="store", dest="mode_selected", type=str)
parser.add_argument('-t', action="store", dest="SNP_info", type=str)
parser.add_argument('-o', action="store", dest="allele_count_output", type=str)
parser.add_argument('-c', action="store", dest="pileup_read_mismatch_threshold", type=float)
parser.add_argument('-d', action="store", dest="min_depth_cov", type=int)
parser.add_argument('-g', action="store", dest="call_haplogroups", type=str)


#test parameters
try:
    options=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

pileup_input=options.pileup_input
mode_selected=options.mode_selected
SNP_info=options.SNP_info
allele_count_output=options.allele_count_output
pileup_read_mismatch_threshold=options.pileup_read_mismatch_threshold

min_depth_cov=options.min_depth_cov

if options.min_depth_cov is None:
    min_depth_cov=0

#if the user specifies a threshold below 0.5, set threshold to 0.5 and print that to screen
if pileup_read_mismatch_threshold<0.5:
    pileup_read_mismatch_threshold=0.5
    print("Using a pileup read mismatch threshold =", str(pileup_read_mismatch_threshold))

# print(pileup_input, mode_selected, allele_count_output, pileup_read_mismatch_threshold)

print ("\nProcessing ", pileup_input)


# pileup_read_mismatch_threshold and number of mismatches
# 0.5 mismatch 11
# 0.6 mismatch 11
# 0.7 mismatch 18
# 0.8 mismatch 26
# 0.9 mismatch 132
# 1 mismatch 406

# mis0.5:7168798 C T 2 1 0
# mis0.7:7168798 C T 2 1 -9
# mis0.9:7168798 C T 2 1 -9

# mis0.5:7204770 G A 3 1 0
# mis0.7:7204770 G A 3 1 0
# mis0.9:7204770 G A 3 1 -9



#################################################################
def remove_deamin(A_count, T_count,Derived_match, mut, mode):
    """
    remove cases where mismateches can be reliably excluded.

    mode 'default' - if mutation is C->T or G->A, then if T_count or A_count==1 (and allele count is = 1),
                          then set A_count=0 or T_count=0.
                          if allele count more than 1, then consider it real,  i.e.
                          T_count or A_count>1, then consider this as a real nutation
    mode 'no-filter' - leave C->T or G->A mutations as they are
    """
    if mut!='C->T' and mut!='G->A':
        if Derived_match=='T':
            if T_count>0:
                T_count=0
                return((T_count))
            else:
                return((T_count))
        elif Derived_match=='A':
            if A_count>0:
                A_count=0
                return((A_count))
            else:
                return((A_count))

    if mut=='C->T' or mut=='G->A':
        if mode=='default':
            if mut=='C->T' and Derived_match=='T':
                if T_count==1:
                    T_count=0
                    return((T_count))
                else:
                    return((T_count))
            elif mut=='G->A' and Derived_match=='A':
                if A_count==1:
                    A_count=0
                    return((A_count))
                else:
                    return((A_count))

        elif mode=='no-filter':
            if mut=='C->T' and Derived_match=='T':
                return(T_count)
            if mut=='G->A' and Derived_match=='A':
                return((A_count))

        # elif mode=='transversions':
        #     if mut=='C->T' and Derived_match=='T':
        #         T_count=0
        #         return(T_count)
        #     if mut=='G->A' and Derived_match=='A':
        #         A_count=0
        #         return(A_count)


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



#this is a good case for deamination
# if REF_allele=='C' and T_count>0:
#         print(variant, A_count,T_count,C_count,G_count)
# ['Y', '16611844', 'C', '2', 't.', 'AF'] 0 1 1 0


#################################################################



#read in haplogroup determination results
pileup=[]
with open(pileup_input, 'r') as in_pileup:
    for line in in_pileup.readlines():
        line=line.strip()
        col=line.split('\t')
        pileup.append(col)
    in_pileup.close()




#################################################################



base_calls = get_bases(pileup)


pos_dict = {}
with open(SNP_info) as f:
    for line in f:
        cols=line.split()
        pos, description = cols[3],line.split()
        pos_dict[pos] = description



res_derived=[]
res_ancestral=[]
res_mismatch=[]
res_removed=[]

output_allele_status=[]

for entry in base_calls:
    CHR = entry[0]
    POS = entry[1]
    REF_allele = entry[2]
    A_count = entry[3]
    T_count = entry[4]
    C_count = entry[5]
    G_count = entry[6]

    if POS in pos_dict:
        chr_match = pos_dict[POS][0]
        marker_match = pos_dict[POS][1]
        hg_match = pos_dict[POS][2]
        pos_match = pos_dict[POS][3]
        mut_match = pos_dict[POS][4]
        Ancestral_match = pos_dict[POS][5].upper()
        Derived_match = pos_dict[POS][6].upper()
        if CHR == chr_match or CHR.replace('chr','') == chr_match:
            if POS == pos_match:
                if mode_selected=='transversions':
                    if mut_match=='C->T' or mut_match=='T->C':
                        T_count=0
                        C_count=0
                    elif mut_match=='G->A' or mut_match=='A->G':
                        A_count=0
                        G_count=0
                else:
                    if mut_match=='C->T' and Derived_match=='T':
                        T_count = remove_deamin(A_count, T_count,Derived_match, mut_match, mode_selected)
                    elif mut_match=='G->A' and Derived_match=='A':
                        A_count = remove_deamin(A_count, T_count,Derived_match, mut_match, mode_selected)
                if eval(Ancestral_match + '_count')>0 and eval(Derived_match + '_count')>0:
                    obs_mismatchAnc = eval(Ancestral_match + '_count')/(eval(Ancestral_match + '_count')+eval(Derived_match + '_count'))
                    obs_mismatchDer = eval(Derived_match + '_count')/(eval(Ancestral_match + '_count')+eval(Derived_match + '_count'))
                    if obs_mismatchAnc!=obs_mismatchDer:
                        if obs_mismatchAnc>obs_mismatchDer and obs_mismatchAnc>=pileup_read_mismatch_threshold:
                            output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 0])
                            res_ancestral.append(hg_match)
                        elif obs_mismatchAnc<obs_mismatchDer and obs_mismatchDer>=pileup_read_mismatch_threshold:
                            output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 1])
                            res_derived.append(hg_match)
                        else:
                            output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), -9])
                            res_mismatch.append(hg_match)
                    else:
                        res_mismatch.append(hg_match)
                        output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), -9])
                elif eval(Derived_match + '_count')>0 and eval(Ancestral_match + '_count')==0:
                    res_derived.append(hg_match)
                    output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 1])
                elif eval(Ancestral_match + '_count')>0 and eval(Derived_match + '_count')==0:
                    res_ancestral.append(hg_match)
                    output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 0])
                elif eval(Ancestral_match + '_count')==0 and eval(Derived_match + '_count')==0:
                    res_removed.append(hg_match)
                    output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), '-9'])
    else:
        print('not in dict')

print("\tread: " + str(len(pileup)) + " calls")
print("\tkept: " + str(len(res_derived)+len(res_ancestral)) + " calls for placement in tree ("+ str(len(res_derived)) +
    " ALTs / "+ str(len(res_ancestral)) + " REFs)")
print("\tset to missing: " + str(len(res_removed)+len(res_mismatch)) + " calls ("+ str(len(res_removed)) +
    " with no data / "+ str(len(res_mismatch)) + " mismatches)")

# print("excluded ", str(len(res_derived)+len(res_ancestral)), "calls ("+ str(len(res_derived)) +
#     " ALTs / "+ str(len(res_ancestral)) + " REFs)")

# print('derived ' + str(len(res_derived)))
# print('ancestral ' + str(len(res_ancestral)))
# print('mismatch ' + str(len(res_mismatch)))
# print('removed ' + str(len(res_removed)))
# print('\n')


#################################################################
with open(allele_count_output, 'w') as out_allele:
    for allele_out in output_allele_status:
        line_to_write=''
        for elem in allele_out:
            line_to_write= line_to_write + str(elem) + ' '
        out_allele.writelines(line_to_write)
        out_allele.writelines('\n')
    out_allele.close()
print("\twritten to " + allele_count_output + '\n\n')


#################################################################
# call_haplogroups routine

if not options.call_haplogroups == 'none':
    print('Calling haplogroups:')
    snps_hg_file=options.call_haplogroups


    hg_pileup=[]
    with open(pileup_input.replace('.pileup','.calls_hgs.pileup'), 'r') as in_hgs_pileup:
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
        print('\tinvalid markers excluded from input haplogroup list: ' + ' ' +  str(excluded_lines))
        print('\tvalid markers for haplogroup calling: ' + ' ' +  str(len(hg_dict)))

    # line of the form ['CHR', 'POS', 'REF','A','T','C','G'], need to keep dict
    base_dict={}
    base_dict['A']=3
    base_dict['T']=4
    base_dict['C']=5
    base_dict['G']=6


    anc_snps=0
    der_snps=0
    NA_snps=0
    with open(allele_count_output.replace('.intree.txt','.hgs.txt'),'w') as outfile:
        hdr = '\t'.join(['POS', 'marker','haplogroup', 'mutation','ANC', 'DER', 'ANC_count', 'DER_count', 'status'])
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
                if ANC_count==0 and DER_count==0:
                    status='NA(nodatapassfilt)'
                    NA_snps=NA_snps+1
                elif ANC_count>0 and DER_count>0:
                    obs_mismatchAnc = ANC_count/ANC_count+DER_count
                    obs_mismatchDer = DER_count/ANC_count+DER_count
                    if obs_mismatchAnc>obs_mismatchDer and obs_mismatchAnc>=pileup_read_mismatch_threshold:
                        status='NA(mismatch)'
                        NA_snps=NA_snps+1
                    elif obs_mismatchDer>obs_mismatchAnc and obs_mismatchDer>=pileup_read_mismatch_threshold:
                        status='NA(mismatch)'
                        NA_snps=NA_snps+1
                else:
                    if (ANC_count>0):
                        status='Anc'
                        anc_snps=anc_snps+1
                    elif (DER_count>0):
                        status='Der'
                        der_snps=der_snps+1

                new_line=[pos, marker, hg, ANC+'->'+DER, ANC, DER,str(ANC_count), str(DER_count), status]
                outfile.writelines('\t'.join(new_line))
                outfile.writelines('\n')
        outfile.close()

        print('\nancestral state: ' + str(anc_snps))
        print('\tderived state: ' + str(der_snps))
        print('\texcluded: ' + str(NA_snps))
        print("\twritten to " + allele_count_output.replace('.intree.txt','.hgs.txt') + '\n\n')

 


    # with open('test_out.txt','w') as outfile:
    #     hdr = ' '.join(['CHR', 'POS', 'REF','A','T','C','G'])
    #     outfile.writelines(hdr)
    #     outfile.writelines('\n')
    #     for line in get_bases(pileup):
    #         new_line=''
    #         for el in line:
    #             if new_line=='':
    #                 new_line = new_line + str(el)
    #             else:
    #                 new_line = new_line + ' ' + str(el)
    #         outfile.writelines(new_line)
    #         outfile.writelines('\n')
    #     outfile.close()


# expect 5 columns
# format:
# CTS7905 NO1     17732587        C       A
# M2330   NO1     17732587        C       A




