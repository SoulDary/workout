from collections import defaultdict


f = input("enter path to fasta: ").strip()
chrs_names = []
new_seqs = []
chrs = []
with open(f) as file:
    while True:
        line = file.readline().strip()
        if line:
            if line.startswith(">"):
                chrs.append("")
                chrs_names.append(line)
            else:
                chrs[-1] += line
        else:
            break

def reverse_kmer(kmer: str):
    complemetary = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    reverse_kmer = ""
    for nucl in reversed(kmer):
        reverse_kmer += complemetary[nucl]
    return reverse_kmer


i = 0
name_i = -1

for seq in chrs:
    name_i += 1
    stikies = defaultdict(list)
    last_site = 0
    i = 0
    while i < len(seq) - 5:
        kmer = seq[i:i+6]
        rev_kmer = reverse_kmer(kmer)
        if kmer == reverse_kmer(kmer):
            first_seq = [seq[last_site:i+1], kmer[1:5]]
            last_site = i+1
            i += 5
            break
        else:
            i += 1

    while i < len(seq) - 5:
        kmer = seq[i:i+6]
        if kmer == reverse_kmer(kmer):
            subseq = seq[last_site:i+1]
            # left_end: forward seq, right_end
            stikies[subseq[:4]].append([subseq, kmer[1:5]])
            last_site = i+1
            i += 5
        else:
            i += 1

    last_seq = [seq[last_site:], kmer[1:5]]


    complemetary = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    new_seq = first_seq[0]


    left_seq = first_seq
    while len(new_seq) < len(seq)-len(last_seq[0]):
        next_seq = None
        nuc_seq = left_seq[0]
        stiky_from_left_seq = left_seq[1]
        for i, subseq in enumerate(stikies[stiky_from_left_seq]):
            if complemetary[subseq[0][4]] != nuc_seq[-1] and subseq[1] in stikies:
                next_seq = stikies[stiky_from_left_seq].pop(i)
                break
        if not next_seq:
            for i, subseq in enumerate(stikies[stiky_from_left_seq]):
                if subseq[1] in stikies:
                    next_seq = stikies[stiky_from_left_seq].pop(i)
                    break
                elif subseq[1] == last_seq[1]:
                    next_seq = subseq[0] + last_seq[0]
                    break
                else:
                    # print("I can't fold it")
                    break
        if len(new_seq) < len(seq)-len(last_seq[0]):
            left_seq = subseq
            new_seq += left_seq[0]
        else:
            break
    #print(new_seq)
    new_seqs.append((new_seq, chrs_names[name_i]))


with open("coronniy.fasta", "w") as output:
    for line in new_seqs:
        output.write(line[1] + "\n")
        output.write(line[0] + "\n")
