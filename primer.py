import math

class Primer:
    def __init__(self, target, start, end):
        self.start = start                  #index of primer start in target seq
        self.end = end                      #index of primer end in target seq
        self.seq = target.seq[start:end]    #primer dna sequence
        self.target = target                #target dna as Target object
        self.len = end-start                #length of primer
    
    def __repr__(self):
        return "5'-"+self.seq+"-3'"

    ## Returns the reverse complement of the primer for self annealing purposes
    def reverse_complement(self):
        complement = {"a":"t", "c":"g", "g":"c", "t":"a"}
        seq = self.seq[::-1]
        rseq = ""
        for letter in seq:
            rseq += complement[letter]
        return rseq
            
    ## Calculate melting temperature of the primer, given salt concentration set to default 50 nM
    def mt(self, salt=0.05):
        #Enthalpy and entropy values for pairs of nucleotides from nearest neighbour thermodynamics model
        nn = {"aa":(9.1,24),"tt":(9.1,24),"at":(8.6,23.9),"ta":(6.0,16.9),\
              "ca":(5.8,12.9),"tg":(5.8,12.9),"gt":(6.5, 17.3),"ac":(6.5, 17.3),\
              "ct":(7.8, 20.8),"ag":(7.8, 20.8),"ga":(5.6, 13.5),"tc":(5.6, 13.5),\
              "cg":(11.9, 27.8),"gc":(11.1, 26.7),"gg":(11.0, 26.6),"cc":(11.0, 26.6)}
        Rlnc4 = -36.1585
        tot = -273.15+(16.6*math.log(salt, 10))
        dH = 0
        dS = 0
        for i in range(self.len-1):
            h, s = nn[self.seq[i:i+2]]
            dH += h
            dS += s
        return ((-dH*1000)/(-dS+Rlnc4))+tot

    ## Calculate gc content of primer as percentage       
    def gc_content(self):
        count = 0
        for letter in self.seq:
            if letter in "gc":
                count += 1
        return count*100/self.len

    ## Returns boolean as to whether GC clamp at 3' end of primer exists
    def gc_clamp(self):
        for pattern in ["cc","cg","gc","gg"]:
            if pattern in self.seq[-5:]:
                return True
        return False

    ## Returns boolean as to whether primer contains a homooligomer of length 5 or greater
    def homooligomer(self):
        for i in range(self.len-5):
            if self.seq[i:i+5] in ["aaaaa","ccccc","ggggg","ttttt"]:
                return True
        return False

    ## Returns boolean as to whether primer contains > 3 x dinucleotide repeat
    def dinucleotide_repeat(self):
        for i in range(self.len-6):
            if self.seq[i:i+2] == self.seq[i+2:i+4] and\
               self.seq[i+2:i+4] == self.seq[i+4:i+6]:
                return True
        return False

    ## Returns boolean as to whether primer is specific based on following criteria:
    ## no more than 3 G/C in last 5 bases or primer and the last 7 bases of primer
    ## with up to 1 mutation doesn't occur anywhere else in target sequence
    def specificity(self):
        gc_count = 0
        gc_count += self.seq[-5:].count("g")
        gc_count += self.seq[-5:].count("c")
        if gc_count > 3:
            return False

        to_match = []
        seq3 = self.seq[-7:]
        for i in range(7):
            for n in ['a','c','g','t']:
                s = seq3[:i] + n + seq3[i+1:]
                to_match.append(s)

        for s in to_match:
            if self.target.seq.count(s) > 1:
                return False
        return True
            

    ## Takes a primer and determines self-alignments with largest
    ## self-annealing and self-end-annealing scores.
    ## Returns those scores along with the alignments
    def self_annealing(self):
        padded_seq = " "*(self.len-1)+self.seq+" "*(self.len-1)
        rc_seq = self.reverse_complement()
        max_sa = 0
        max_sea = 0
        sa_alignments = []
        sea_alignments = []

        for i in range(2*self.len-1): ##for each alignment
            score = 0
            matches = ""
            for j in range(self.len): ##for each base in primer
                if padded_seq[i+j] == rc_seq[j]:
                    matches += "|"
                    if rc_seq[j] in "at":
                        score += 2
                    else:
                        score += 4
                else:
                    matches += " "

            if i > self.len:
                alignment = "5'-"+self.seq+"-3'"+"\n"\
                            + " "*(i-self.len+4)+matches + "\n"\
                            +" "*(i-self.len+1)+"3'-"+self.seq[::-1]+"-5'\n"
            else:
                alignment = " "*(self.len-i-1)+"5'-"+self.seq+"-3'"+"\n"\
                            +"   "+matches + "\n"\
                            +"3'-"+self.seq[::-1]+"-5'\n"
                
            if score > max_sa:
                max_sa = score
                sa_alignments = [alignment]             
            elif score == max_sa:
                sa_alignments.append(alignment)

            if i >= self.len:
                if score > max_sea:
                    max_sea = score
                    sea_alignments = [alignment]
                elif score == max_sea:
                    sea_alignments.append(alignment)
            
        return (max_sa, "\n".join(sa_alignments), max_sea, "\n".join(sea_alignments))


    ## Takes a forward and reverse primer pair and determines alignments with largest
    ## pair-annealing and pair-end-annealing scores.
    ## Returns those scores along with the alignments
    def pair_annealing(self, other):
        padded_seq = " "*(other.len-1)+self.seq+" "*(other.len-1)
        rc_seq = other.reverse_complement()
        max_pa = 0
        max_pea = 0
        pa_alignments = []
        pea_alignments = []

        for i in range(self.len+other.len-1): ##for each alignment
            score = 0
            matches = ""
            for j in range(other.len): ##for each base in primer
                if padded_seq[i+j] == rc_seq[j]:
                    matches += "|"
                    if rc_seq[j] in "at":
                        score += 2
                    else:
                        score += 4
                else:
                    matches += " "

            if i > self.len:
                alignment = "5'-"+self.seq+"-3'"+"\n"\
                            + " "*(i-self.len+4)+matches + "\n"\
                            +" "*(i-self.len+1)+"3'-"+other.seq[::-1]+"-5'\n"
            else:
                alignment = " "*(other.len-i-1)+"5'-"+self.seq+"-3'"+"\n"\
                            +"   "+matches + "\n"\
                            +"3'-"+other.seq[::-1]+"-5'\n"
            
            if score > max_pa:
                max_pa = score
                pa_alignments = [alignment]             
            elif score == max_pa:
                pa_alignments.append(alignment)

            if i >= self.len:
                if score > max_pea:
                    max_pea = score
                    pea_alignments = [alignment]
                elif score == max_pea:
                    pea_alignments.append(alignment)
        return (max_pa, "\n".join(pa_alignments), max_pea, "\n".join(pea_alignments))
    

    ## Returns the following primer parameters: mt, gc content, gc_clamp
    ## homooligomerisation, dinucleotide repeat, specificity
    def parameters(self):
        return (self.mt(), self.gc_content(), self.gc_clamp(), self.homooligomer(),\
                self.dinucleotide_repeat(), self.specificity())

    ## Takes a forward and reverse primer pair and returns scoring vector:
    ## length, gc_content, mt, sa and sea of both, and pa and pea scores
    def scoring_vector(self, other):
        sa_p, sa_pa, sea_p, sea_pa = self.self_annealing()
        sa_q, sa_qa, sea_q, sea_qa = other.self_annealing()
        pa, pa_a, pea, pea_a = self.pair_annealing(other)
        return (self.len, other.len, self.gc_content(), other.gc_content(),\
                self.mt(), other.mt(), sa_p, sa_q, sea_p, sea_q, pa, pea)



class Probe:
    def __init__(self, target):
        pass


class Target:
    def __init__(self, seq, start, end, is_forward):
        self.seq = seq                  #sequence of entire dna fragment
        self.start = start              #index of start of desired target seq
        self.end = end                  #index of end of desired target seq
        self.len = end-start            #legnth of target sequence
        self.is_forward = is_forward    #boolean as to whether strand is the forward strand

    def __repr__(self):
        fseq = self.seq
        fseq = fseq[:self.start]+fseq[self.start:self.end].upper()+fseq[self.end:]
        
        rstrand = self.reverse_complement()
        rseq = rstrand.seq
        rseq = rseq[:rstrand.start]+rseq[rstrand.start:rstrand.end].upper()+rseq[rstrand.end:]
        rseq = rseq[::-1]
        
        chunks = []
        count = 0
        while len(fseq) > 100:
            chunks.append(str(count*100+1)+" 5'-"+fseq[:100]+"-3' "+str(count*100+100)+"\n"+
                          " "*len(str(count*100+1))+" 3'-"+rseq[:100]+"-5'")
            fseq = fseq[100:]
            rseq = rseq[100:]
            count += 1
        if len(fseq) != 0:
            chunks.append(str(count*100+1)+" 5'-"+fseq+"-3' "+str(len(self.seq))+"\n"+
                          " "*len(str(count*100+1))+" 3'-"+rseq+"-5'")
        return "\n\n".join(chunks)+"\n"
        
    def reverse_complement(self):
        complement = {"a":"t", "c":"g", "g":"c", "t":"a"}
        seq = self.seq[::-1]
        rseq = ""
        for letter in seq:
            rseq += complement[letter]
        start = len(self.seq) - self.end
        end = len(self.seq) - self.start
        is_forward = not self.is_forward
        rc = Target(rseq, start, end, is_forward)
        return rc

    def gc_content(self):
        count = 0
        for letter in self.seq[self.start:self.end]:
            if letter in "gc":
                count += 1
        return count*100/(self.end-self.start)
               



