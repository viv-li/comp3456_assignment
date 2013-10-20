#!/usr/bin/env python3.3
from primer_search import *
import sys

## Get sequence and target from stdin or prompt for it
args = sys.argv
if len(args) > 3:
    sequence = clean_dna("".join(open(args[1],"rU").read().split("\n")[1:]))
    start = int(args[2])
    end = int(args[3])
else:
    file = input("Enter FASTA file path or sequence: ")
    if is_dna(file):
        sequence = clean_dna(file)
    else:
        sequence = clean_dna("".join(open(file,"rU").read().split("\n")[1:]))    

    start = int(input("Enter start of target sequence: "))
    end = int(input("Enter end of target sequence: "))
    if start >= end or end > len(sequence):
        raise Exception("Incorrect start or end given")

fstrand = Target(sequence, start, end, True)
rstrand = fstrand.reverse_complement()


## Get mt, gc content, gc clamp parameters or use defaults

mt = (55,65)
gc_content = (fstrand.gc_content()-10,fstrand.gc_content()+15)
gc_clamp = None

extra = input("Set optional MT, GC content and GC clamp parameters? (Y/N): ")
if extra:
    if extra[0].upper() == "Y":
        user_mt = input("Enter MT range (hit enter for default 57-63°C): ")
        if user_mt:
            mt = tuple(map(int, map(strip, user_mt.split("-"))))
        
        user_gc_content = input("Enter GC content (hit enter for default): ")
        if user_gc_content:
            gc_content = tuple(map(int, map(strip, user_gc_content.split("-"))))
            
        user_gc_clamp = input("Include GC clamp? (Y/N) (hit enter for both): ")
        if user_gc_clamp:
            if user_gc_clamp[0].upper() == "Y":
                gc_clamp = True
            elif user_gc_clamp[0].upper() == "N":
                gc_clamp = False


## Get output location
output = input("Enter location to save results or hit enter to print to stdout: ")


## Find possible primers
user_parameters = (mt, gc_content, gc_clamp)
search_range = 50
while search_range < 300:
    fprimers = find_primers(fstrand, user_parameters, search_range)
    if fprimers:
        break
    search_range += 50
    
search_range = 50
while search_range < 300:
    rprimers = find_primers(rstrand, user_parameters, search_range)
    if rprimers:
        break
    search_range += 50
results = score_pairs(fprimers, rprimers)

## Output results
if not output:
    output = sys.stdout
else:
    output = open(output, "w")

if not results:
    print("No primer pairs found")
for i in range(len(results)):
    result = results[i]
    fprimer = fprimers[result[1]]
    rprimer = rprimers[result[2]]
    sv = fprimer.scoring_vector(rprimer)
    Ta = min(sv[4],sv[5])-5

    print("{:d}. Forward primer: {:s} start = {:d}, length = {:d}bp, %GC = {:.2f}%, Tm = {:.2f}°C"\
          .format(i+1,fprimer,fprimer.start,sv[0],sv[2],sv[4]), file=output)
    print("   Reverse primer: {:s} start = {:d}, length = {:d}bp, %GC = {:.2f}%, Tm = {:.2f}°C"\
          .format(rprimer,len(fstrand.seq)-rprimer.start,sv[1],sv[3],sv[5]), file=output)
    print("   Product size = {:d}bp, Optimal annealing temperature = {:.2f}°C\n"\
          .format(len(fstrand.seq)-rprimer.start-fprimer.start, Ta),file=output)

#print(fstrand, file=output)

if output != sys.stdout:
    output.close()
