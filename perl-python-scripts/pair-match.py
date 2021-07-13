#!/usr/bin/env python
from Bio import AlignIO
import argparse
import sys

parser = argparse.ArgumentParser(description='Pairwise seq diff in sliding windows. Authors: Roman, Edagar, Weigang')
parser.add_argument('fasta', type=str, help='Input file (FASTA alignment)')
parser.add_argument('-s', '--win_size', type=int, default=6,
                                         help='Window size (default = 6)' )
parser.add_argument('-p', '--win_step', type=int, default=1,
                                        help='Window step (dfault = 1')
parser.add_argument('-c', '--cut_off', type=float, default=0,
                                        help='cutoff for minimum avg diff in a window (default = 0)')
parser.add_argument('-w', '--wide', action = 'store_true',
                                        help='print in wide format')

args = parser.parse_args()
align = AlignIO.read(args.fasta, "fasta")
win_size = args.win_size
win_step = args.win_step
aln_len = align.get_alignment_length()
seq_num = len(align)
cut = args.cut_off

pol = {"A":7,"C":4.8,"D":13,"E":12.5,"F":5,"G":7.9,"H":8.4,"I":4.9,"K":10.1,"L":4.9,"M":5.3,"N":10,"P":6.6,"Q":8.6,"R":9.1,"S":7.5,"T":6.6,"V":5.6,"W":5.2,"Y":5.4}
hydro = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
vol = {"A":31,"C":55,"D":54,"E":8.3,"F":132,"G":3,"H":96,"I":111,"K":119,"L":111,"M":105,"N":56,"P":32.5,"Q":85,"R":124,"S":32,"T":61,"V":84,"W":170,"Y":136}
iso = {"A":6,"C":5.07,"D":2.77,"E":3.22,"F":5.48,"G":5.97,"H":7.59,"I":6.02,"K":9.74,"L":5.98,"M":5.74,"N":5.41,"P":6.3,"Q":5.65,"R":10.76,"S":5.68,"T":6.16,"V":5.96,"W":5.89,"Y":5.66}

slices = {}
slice_selected = []

def pair_diff(win, aln):
    matches = []
    for i in range(seq_num):
        seq_i = aln[i,:] # all columns
        for j in range(i+1, seq_num):
            seq_j = aln[j,:]
            match = 0
            diff = 0
            valid = 0
            diffPol = 0
            diffHyd = 0
            diffVol = 0
            diffIso = 0
            for k in range(aln.get_alignment_length()):
                if seq_i[k] == "-" or seq_j[k] == "-":
                    continue
                else:
                    valid += 1
                    if seq_i[k] == seq_j[k]:
                        match += 1
                    else:
                        diff += 1
                        diffPol += abs(pol[seq_i[k]] - pol[seq_j[k]])
                        diffHyd += abs(hydro[seq_i[k]] - hydro[seq_j[k]])
                        diffVol += abs(vol[seq_i[k]] - vol[seq_j[k]])
                        diffIso += abs(iso[seq_i[k]] - iso[seq_j[k]])

            if valid > 0:        
                matches.append({
                    'win': win,
                    "id1": seq_i.id,
                    'seq1': seq_i.seq,
                    "id2": seq_j.id,
                    'seq2': seq_j.seq,
                    'valid': valid,
                    'match': match,
                    'diff_raw': round(1-match/valid,4),
                    'diff_pol': round(diffPol/valid,4),
                    'diff_hydro': round(diffHyd/valid,4),
                    'diff_vol': round(diffVol/valid,4),
                    'diff_iso': round(diffIso/valid,4)                
                })
    slices[win] = matches

def rm_window_by_avg_diff():
    for win in slices:
        sumDiff = 0
        ct = 0
        for s in slices[win]:
            ct += 1
            sumDiff += s['diff_raw']
        if sumDiff/ct >= cut:
            slice_selected.append(slices[win])
        
def slide_windows():
    for i in range(aln_len - win_size + 1):
        aln_slice = align[:, i:i+win_size] # all rows
        pair_diff(i, aln_slice)   
        
slide_windows()
rm_window_by_avg_diff()

if args.wide:
    print("to be implemented, but not be necessary")
else:
    print("\t".join(["window", "id1", "seq1", "id2", "seq2", "valid", "match", 'diff_raw', 'diff_pol', 'diff_hydro', 'diff_vol', "diff_iso"]))
    for matches in slice_selected:
        for match in matches:
            print("\t".join(map(str, [match['win'],
                                  match['id1'],
                                  match['seq1'],
                                  match['id2'],
                                  match['seq2'],
                                  match['valid'],
                                  match['match'],
                                  match['diff_raw'],
                                  match['diff_pol'],
                                  match['diff_hydro'],
                                  match['diff_vol'],
                                  match['diff_iso']
            ])))

#print(len(slice_selected))
#for s in slice_selected:

#

sys.exit()
