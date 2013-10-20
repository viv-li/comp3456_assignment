#!/usr/bin/env python3.3
import primer_search
import cgi

fields = cgi.FieldStorage()
if fields.has_key('sequence')
    sequence = fields['sequence'].value

sequence = clean_dna("".join(open("hiv.fasta","rU").read().split("\n")[1:]))
fstrand = Target(sequence, 5771, 8341, True)
rstrand = fstrand.reverse_complement()

##mt, gc content, gc clamp defaults
user_parameters = ((57,63), (fstrand.gc_content()-5,fstrand.gc_content()+10), None)
fprimers = find_primers(fstrand, user_parameters)
rprimers = find_primers(rstrand, user_parameters)

print('Content-Type: text/html\n')
print('''<!DOCTYPE html>
<html>
  <head>
    <title>Hello %s!!</title>
  </head>
  <body>''')

count = 1
for result in score_pairs(fprimers, rprimers):
    fprimer = fprimers[result[1]]
    fself = fprimer.self_annealing()
    rprimer = rprimers[result[2]]
    rself = rprimer.self_annealing()
    pair = fprimer.pair_annealing(rprimer)
    sv = fprimer.scoring_vector(rprimer)
    Ta = min(sv[4],sv[5])-5

    print("{:d}. Forward primer: {:s} start = {:d}, length = {:d}bp, %GC = {:f}, Tm = {:.2f}°C".format(count,fprimer,fprimer.start,sv[0],sv[2],sv[4]))
    print("   Reverse primer: {:s} start = {:d}, length = {:d}bp, %GC = {:f}, Tm = {:.2f}°C".format(rprimer,len(fstrand.seq)-rprimer.start,sv[1],sv[3],sv[5]))
    print("   Product size = {:d}bp, Ta = {:.2f}°C".format(len(fstrand.seq)-rprimer.start-fprimer.start, Ta))
##    print("\n".join(map(str,pair)))
    print()
    count += 1

print('''</body>
</html>''')