#! /usr/env python3

from joblib import Parallel, delayed
import sys, gzip, subprocess
from tempfile import TemporaryDirectory
import traceback


def findOp (cigarStr, length, op = "I"):
    """Return index of most similar op"""
    from cigar import Cigar
    c = Cigar(cigarStr)
    score = length
    n = -1
    result = n
    for span, opName in c.items():
        n += 1
        if opName != op:
            continue
        diff = abs(span - length)
        if diff < score:
            score = diff
            result = n
    return result

def assemble (fasta, tmpDir, id, insSeq, length):
    from Bio.pairwise2 import align, format_alignment
    asmCommand = ["wtdbg2", "-i", fasta.name, "-t", "1", "-L", "1000",
        "-o", tmpDir + "/" + id, "-x", "ont"
    ]
    try:
        asmProcess = subprocess.run(asmCommand, capture_output = True, 
            check = True,
            timeout = 6 * 60)
    except:
        return (insSeq, "LENGTH")
    consCommand = ["wtdbg-cns", "-t", "1", "-i", 
        tmpDir + "/" + id + ".ctg.lay.gz", 
    ]
    try:
        cons = subprocess.run(consCommand, capture_output = True, 
            check = True,
            timeout = 6 * 60)
    except:
        return (insSeq, "LENGTH")

    asmResult = cons.stdout.decode()
    asm = "".join(asmResult.rstrip().split("\n")[1:])
    if asmResult == "":
        return (insSeq, "LENGTH")
    aln = align.localxs(asm, insSeq, -10, -0.5, one_alignment_only = True)[0]
    result = asm[aln.start:aln.end]
    return (result, "ASM")

vcfFile = sys.argv[1]
bamFile = sys.argv[2]
outputFile = sys.argv[3]
outputDir = "."

out = []

def processSv (line, parentTmp, bamFile):
    from math import floor, ceil
    from pysam import AlignmentFile
    from tempfile import TemporaryDirectory
    import re
    sv = line.split("\t")
    chrom = sv[0]
    pos = int(sv[1])
    idSV = sv[2]
    info = sv[7]
    with TemporaryDirectory(prefix = idSV, dir = parentTmp) as individualTmp, \
        AlignmentFile(bamFile) as bam,\
        open(individualTmp + "/fasta.fasta", "w") as fasta:
        support = re.search(r"(?:RE|SUPPORT)=([^;\t]+)", info).groups()[0]
        if int(support) < 5:
            sv[7] += ";ALTMETHOD=CALLER"
            return("\t".join(sv))
        readNames = re.search(r"(?:RNAMES|READS)=([^;\t]+)", info).groups()[0]
        if readNames == "NULL":
            sv[7] += ";ALTMETHOD=CALLER"
        else:
            try:
                readNames = readNames.split(",")
                svLen = abs(int(re.search(r"SVLEN=([^;\t]+)", info).groups()[0]))
                reads = bam.fetch(chrom, pos - 100, pos + 100)
                sequences = []
                positions = []
                usefulReads = 0
                for r in reads:
                    if r.query_name not in readNames or r.mapping_quality < 20 or not r.query_sequence:
                        continue
                    bestOp = "" # the tuple ("OP", length) best matching the insertions
                    bestDiff = svLen
                    insIndex = 0
                    for i in range(len(r.cigartuples)):
                        opCode, opLength = r.cigartuples[i]
                        if opCode != 1:
                            continue
                        difference = abs(opLength - svLen)
                        if difference < bestDiff:
                            insIndex = i
                            bestDiff = difference
                            bestOp = (opCode, opLength)
                    if bestOp[1] not in range(floor(svLen * 0.85), ceil(svLen * 1.15)):
                        continue
                    startInRead = 0
                    referenceShift = 0
                    for op, length in r.cigartuples[0:insIndex]: # the insertion itself is not counted
                        if op in [0, 7, 8]: # M, =, X
                            startInRead += length
                            referenceShift += length
                        elif op in [4, 1]: #S, I
                            startInRead += length
                        elif op in [2, 3]: # D, N
                            referenceShift += length
                        # H (5) and P (6) dont affect anything
                    endInRead = startInRead + r.cigartuples[insIndex][1]
                    beginning = 0
                    end = r.query_length
                    for i in range(len(r.cigartuples)):
                        tup = r.cigartuples[i]
                        if tup[0] == 4:
                            beginning += tup[1]
                        elif tup[0] == 5:
                            continue
                        else:
                            break
                    i = len(r.cigartuples) - 1
                    while i >= 0:
                        tup = r.cigartuples[i]
                        if tup[0] == 4:
                            end += tup[1]
                            i -= 1
                        elif tup[0] == 5:
                            i -= 1
                            continue
                        else:
                            break
                    seq = r.query_sequence[startInRead:endInRead]
                    positions.append(referenceShift + r.reference_start + 1)
                    sequences.append(seq)
                    if startInRead - 5000 < beginning:
                        compensation = abs(beginning - abs(startInRead - 5000))
                        seq = r.query_sequence[beginning:min(endInRead + 5000 + compensation, end)]
                    elif endInRead + 5000 > end:
                        compensation = endInRead + 5000 - end
                        seq = r.query_sequence[max(startInRead - 5000 - compensation, beginning):end]
                    else:
                        seq = r.query_sequence[max(startInRead - 5000, beginning):min(endInRead + 5000, end)]
                    if len(seq) < 1000:
                        continue
                    usefulReads += 1
                    fasta.write(f">{r.query_name}\n{seq}\n")
                fasta.close()
                if len(sequences) < 5 or usefulReads < 5:
                    sv[7] += ";ALTMETHOD=CALLER"
                    return "\t".join(sv)
                lengthDiff = svLen
                bestIns = ""
                for i in sequences:
                    difference = abs(len(i) - svLen)
                    if difference < lengthDiff:
                        lengthDiff = difference
                        bestIns = i
                        if difference == 0:
                            break
                alt, method = assemble(fasta, 
                    individualTmp, 
                    idSV, bestIns, svLen)
                sv[4] = alt
                sv[7] = re.sub("(RNAMES|READS)=[^;]+", "", sv[7]) + "ALTMETHOD=" + method
            except BaseException as e:
                traceback.print_exc()
                print("Exception for ", "\t".join(sv))
                return "\t".join(sv)
        return("\t".join(sv))

with TemporaryDirectory(prefix="polishINS_") as parentTmp:
    with gzip.open(vcfFile, mode = "rt") as vcf, \
        open(outputFile, "w") as outVcf:

        n = 0
        while 1: # Lets read the header
            line = vcf.readline().rstrip()
            if line.startswith("##"):
                out.append(line)
            elif line.startswith("#C"):
                out.append('##INFO=<ID=ALTMETHOD,Number=1,Type=String,Description="Obtention of ALT">')
                out.append(line)
                break
        outVcf.write("\n".join(out))
        outVcf.write("\n")
        out = []
        lines = [x.rstrip() for x in vcf.readlines() if "SVTYPE=INS" in x]
        print("INS number: ", len(lines))
        out = Parallel(n_jobs = 15)(delayed(processSv)(x, parentTmp, bamFile) for x in lines)

        outVcf.write("\n".join(out))
        outVcf.write("\n")
