from primer import Primer
from primer import Target

def weighted_distance(sc, sc_ideal):
    wd = 0
    ##distances with weights applied
    wd += abs(sc[0] - sc[1])*0.5 ##equal length
    wd += abs(sc[2] - sc[3]) ##equal gc content
    wd += abs(sc[4] - sc[5]) ##equal melting temp
    wd += abs(sc[6] - sc_ideal[6])*0.1 ##forward sa
    wd += abs(sc[7] - sc_ideal[7])*0.1 ##reverse sa
    wd += abs(sc[8] - sc_ideal[8])*0.2 ##forward sea
    wd += abs(sc[9] - sc_ideal[9])*0.2 ##reverse sea
    wd += abs(sc[10] - sc_ideal[10])*0.1 ##pa
    wd += abs(sc[11] - sc_ideal[11])*0.2 ##pea

    return wd


def find_primers(target_dna, user_parameters):
    primers = []
    for i in range(target_dna.start+1):
        for j in range(12, 21):
            primer = Primer(target_dna, i, i+j)
            parameters = primer.self_score()
            
            if user_parameters[0][0] <= parameters[0] <= user_parameters[0][1] and\
               user_parameters[1][0] <= parameters[1] <= user_parameters[1][1]:
                if user_parameters[2] == 1:
                    if parameters[2] ==1:
                        primers.append(primer)
                else:
                    primers.append(primer)
                    
    return primers

def score_pairs(fprimers, rprimers, ideal_sv):
    scores = []
    i = 0
    for f in fprimers:
        j = 0
        for r in rprimers:
            scores.append((weighted_distance(f.scoring_vector(r), ideal_sv), i, j))
            j += 1
        i += 1
    scores.sort()
    return scores[:5]


fstrand = Target("gattaagtttctagtcaccgctcgatatgctagctagactggctacg", 10, 30)
rstrand = fstrand.reverse_complement()
print(fstrand)
print(rstrand)

user_parameters = ((50,60), (40,60), 0)
fprimers = find_primers(fstrand, user_parameters)
rprimers = find_primers(rstrand, user_parameters)
##print(fprimers)
##print(rprimers)

ideal_sv = (18, 18, 50, 50, 55, 55, 0,0,0,0,0,0)

for score in score_pairs(fprimers, rprimers, ideal_sv):
    print(score[0], fprimers[score[1]], rprimers[score[2]], fprimers[score[1]].scoring_vector(rprimers[score[2]]))

##mt, gc content, gc clamp defaults
min_len = 12
max_len = 20

 
