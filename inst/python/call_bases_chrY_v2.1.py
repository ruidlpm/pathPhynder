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
parser.add_argument('-p', action="store", dest="proportion_of_mismatches_tolerated", type=float)


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
if options.proportion_of_mismatches_tolerated:
    proportion_of_mismatches_tolerated=options.proportion_of_mismatches_tolerated
else:
    proportion_of_mismatches_tolerated=0.3
    print(pileup_input, mode_selected, allele_count_output, proportion_of_mismatches_tolerated)







#################################################################
def remove_deamin(A_count, T_count,Derived_match, mut, mode):
    """
    remove cases where mismateches can be reliably excluded.
    mode 'conservative' - if mutation is C->T or G->A, then if T_count or A_count==1 (and allele count is = 1),
                          then set A_count=0 or T_count=0.
                          if allele count more than 1, then consider it real,  i.e.
                          T_count or A_count>1, then consider this as a real nutation
    mode 'relaxed' - leave C->T or G->A mutations as they are
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
        if mode=='conservative':
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

        elif mode=='relaxed':
            if mut=='C->T' and Derived_match=='T':
                return(T_count)
            if mut=='G->A' and Derived_match=='A':
                return((A_count))




#################################################################



#read in haplogroup determination results
pileup=[]
with open(pileup_input, 'r') as in_pileup:
    for line in in_pileup.readlines():
        line=line.strip()
        col=line.split('\t')
        pileup.append(col)
    in_pileup.close()

print("read", len(pileup), "calls")



#################################################################

def chunk_string(s, n):
    return [s[i:i+n] for i in range(len(s)-n+1)]


def get_bases(pileup):
    #first line is header - add hg(if any) and branch(later)
    parsed_pileup=[]
    for variant in pileup:
        CHR = variant[0]
        POS = variant[1]
        #clear '$' and '^*' patterns
        
        if variant[3] =='0':
            variant.append('0')
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
    #     if len(variant[4]) != len(variant[5]):
        #count reference bases
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

with open('test_out.txt','w') as outfile:
    hdr = ' '.join(['CHR', 'POS', 'REF','A','T','C','G'])
    outfile.writelines(hdr)
    outfile.writelines('\n')
    for line in get_bases(pileup):
        new_line=''
        for el in line:
            if new_line=='':
                new_line = new_line + str(el)
            else:
                new_line = new_line + ' ' + str(el)
        outfile.writelines(new_line)
        outfile.writelines('\n')
    outfile.close()


base_calls = get_bases(pileup)
#################################################################

pos_dict = {}
with open(SNP_info) as f:
    for line in f:
        cols=line.split()
        pos, description = cols[3],line.split()
        pos_dict[pos] = description




res_derived=[]
res_ancestral=[]
res_mismatch=[]

output_allele_status=[]

for entry in base_calls:
    CHR = entry[0]
    POS = entry[1]
    REF_allele = entry[2]
    A_count = entry[3]
    T_count = entry[4]
    C_count = entry[5]
    G_count = entry[6]

    entry_match = pos_dict[POS]
    chr_match = pos_dict[POS][0]
    marker_match = pos_dict[POS][1]
    hg_match = pos_dict[POS][2]
    pos_match = pos_dict[POS][3]
    mut_match = pos_dict[POS][4]
    Ancestral_match = pos_dict[POS][5].upper()
    Derived_match = pos_dict[POS][6].upper()
    if CHR == chr_match or CHR.replace('chr','') == chr_match:
        if POS == pos_match:
#            if mut_match!='C->T' and Derived_match=='T':
#                T_count = remove_deamin(A_count, T_count,Derived_match, mut_match, mode_selected)
#            elif mut_match!='G->A' and Derived_match=='A':
#                A_count = remove_deamin(A_count, T_count,Derived_match, mut_match, mode_selected)
            if mut_match=='C->T' and Derived_match=='T':
                T_count = remove_deamin(A_count, T_count,Derived_match, mut_match, mode_selected)
            elif mut_match=='G->A' and Derived_match=='A':
                A_count = remove_deamin(A_count, T_count,Derived_match, mut_match, mode_selected)
            #simple evaluation if derived or ancestral alleles. Add also a 1 if derived, 0 if ancestral, -9 if mismatch
            # note, prob better to add -9 if Anc==Der, and estimate prop Der/Anc for cases where Anc!=Der.
            #for now only -9
            if eval(Ancestral_match + '_count')>0 and eval(Derived_match + '_count')>0:
                obs_mismatchAnc = eval(Ancestral_match + '_count')/(eval(Ancestral_match + '_count')+eval(Derived_match + '_count'))
                obs_mismatchDer = eval(Derived_match + '_count')/(eval(Ancestral_match + '_count')+eval(Derived_match + '_count'))
                if obs_mismatchAnc!=obs_mismatchDer:
                    if obs_mismatchAnc>obs_mismatchDer and obs_mismatchAnc>=(1-proportion_of_mismatches_tolerated):
                        output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 0])
                    elif obs_mismatchAnc<obs_mismatchDer and obs_mismatchDer>=(1-proportion_of_mismatches_tolerated):
                        output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 1])
                    else:
                        output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), -9])
                else:
                    res_mismatch.append(hg_match)
                    output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), -9])
            elif eval(Derived_match + '_count')>0 and eval(Ancestral_match + '_count')==0:
                res_derived.append(hg_match)
                output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 1])
            elif eval(Ancestral_match + '_count')>0 and eval(Derived_match + '_count')==0:
                res_ancestral.append(hg_match)
                output_allele_status.append([POS, Ancestral_match, Derived_match, eval(Ancestral_match + '_count'), eval(Derived_match + '_count'), 0])

#################################################################
with open(allele_count_output, 'w') as out_allele:
    for allele_out in output_allele_status:
        line_to_write=''
        for elem in allele_out:
            line_to_write= line_to_write + str(elem) + ' '
        out_allele.writelines(line_to_write)
        out_allele.writelines('\n')
    out_allele.close()
