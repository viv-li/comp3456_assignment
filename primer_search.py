from primer import Primer
from primer import Target

def clean_dna(seq):
    cleaned = ""
    for letter in seq.lower():
        if letter in "acgt":
            cleaned += letter
    return cleaned

#given the scoring vector 
def weighted_distance(sv):
    wd = 0
    ##distances with weights applied
    wd += abs(sv[0] - sv[1])*0.5 ##equal length
    wd += abs(sv[2] - sv[3]) ##equal gc content
    wd += abs(sv[4] - sv[5]) ##equal melting temp
    wd += sv[6]*0.1 ##forward sa
    wd += sv[7]*0.1 ##reverse sa
    wd += sv[8]*0.2 ##forward sea
    wd += sv[9]*0.2 ##reverse sea
    wd += sv[10]*0.1 ##pa
    wd += sv[11]*0.2 ##pea
    return wd


## given a single strand of DNA, find all acceptable primers
## mt, gc content, gc clamp or no, no homooligomers, no dinucleotide repeats, specificity
def find_primers(target_dna, user_parameters):
    primers = []
    if target_dna.len < 500:
        ideal_len = (16, 18)
        search_range = 50
    elif target_dna.len > 2500:
        ideal_len = (20, 24)
        search_range = 250
    else:
        ideal_len = (18, 20)
        search_range = 100
    search_start = max(0, target_dna.start-search_range)
    
    for i in range(search_start, target_dna.start+1):
        for j in range(ideal_len[0], ideal_len[1]+1):
            primer = Primer(target_dna, i, i+j)
            parameters = primer.parameters()
            
            if user_parameters[0][0] <= parameters[0] <= user_parameters[0][1] and\
               user_parameters[1][0] <= parameters[1] <= user_parameters[1][1] and\
               not parameters[3] and not parameters[4] and parameters[5]:
                if user_parameters[2] == None:
                    primers.append(primer)
                else:
                    if user_parameters[2] == parameters[2]:
                        primers.append(primer)                   
    return primers

## given lists of forward and reverse primers, pair each up and calculate
## scoring vector, return 5 best pairs
def score_pairs(fprimers, rprimers):
    scores = []
    i = 0
    for f in fprimers:
        j = 0
        for r in rprimers:
            scores.append((weighted_distance(f.scoring_vector(r)), i, j))
            j += 1
        i += 1
    scores.sort()
    return scores[:5]
 
