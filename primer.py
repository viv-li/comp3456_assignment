import math

class Primer:
    def __init__(self, target, start=0, end=18):
        self.start = start
        self.end = end
        self.seq = target.seq[start:end]
        self.len = end-start
    
    def __repr__(self):
        return "5'-"+self.seq+"-3'"

    def reverse_complement(self):
        complement = {"a":"t", "c":"g", "g":"c", "t":"a"}
        seq = self.seq[::-1]
        rseq = ""
        for letter in seq:
            rseq += complement[letter]
        return rseq
            
    ## Calculate melting temperature of the primer, given salt concentration set to default 50 nM
    def mt(self, salt=0.05):
        Rlnc4 = -36.1585
        tot = -273.15+(16.6*math.log(salt, 10))
        nn = {"aa":(9.1,24),"tt":(9.1,24),"at":(8.6,23.9),"ta":(6.0,16.9),\
              "ca":(5.8,12.9),"tg":(5.8,12.9),"gt":(6.5, 17.3),"ac":(6.5, 17.3),\
              "ct":(7.8, 20.8),"ag":(7.8, 20.8),"ga":(5.6, 13.5),"tc":(5.6, 13.5),\
              "cg":(11.9, 27.8),"gc":(11.1, 26.7),"gg":(11.0, 26.6),"cc":(11.0, 26.6)}

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
        if self.seq[-1] in "gc" and self.seq[-2] in "gc":
            return True
        return False

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
                            +" "*(i-self.len+1)+"3'-"+self.seq[::-1]+"-5'"
            else:
                alignment = " "*(self.len-i-1)+"5'-"+self.seq+"-3'"+"\n"\
                            +"   "+matches + "\n"\
                            +"3'-"+self.seq[::-1]+"-5'"
                
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

##        for a in sa_alignments:
##            print(a)
##        for a in sea_alignments:
##            print(a)
            
        return (max_sa, max_sea)

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
                            +" "*(i-self.len+1)+"3'-"+other.seq[::-1]+"-5'"
            else:
                alignment = " "*(other.len-i-1)+"5'-"+self.seq+"-3'"+"\n"\
                            +"   "+matches + "\n"\
                            +"3'-"+other.seq[::-1]+"-5'"
            
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

##        for a in pa_alignments:
##            print(a)
##        for a in pa_alignments:
##            print(a)
            
        return (max_pa, max_pea)        

    ## Takes a forward and reverse primer pair and
    ## returns their scoring vector
    def scoring_vector(self, other):
        sa_p, sea_p = self.self_annealing()
        sa_q, sea_q = other.self_annealing()
        pa, pea = self.pair_annealing(other)
        return (self.len, other.len, self.gc_content(), other.gc_content(),\
                self.mt(), other.mt(), sa_p, sa_q, sea_p, sea_q, pa, pea)
        
    def self_score(self):
        return (self.mt(), self.gc_content(), self.gc_clamp())
        

class Target:
    def __init__(self, seq, start, end):
        self.seq = seq
        self.start = start
        self.end = end
        self.len = len(seq)

    def __repr__(self):
        return "5'-" + self.seq[:self.start]+self.seq[self.start:self.end].upper()\
               +self.seq[self.end:] + "-3'"
        
    def reverse_complement(self):
        complement = {"a":"t", "c":"g", "g":"c", "t":"a"}
        seq = self.seq[::-1]
        rseq = ""
        for letter in seq:
            rseq += complement[letter]
        start = self.len - self.end
        end = self.len - self.start
        rc = Target(rseq, start, end)
        return rc

               



