#!/bin/env python

import subprocess
import sys
import argparse
import pathlib


def writeIndexFqs(inpFq, outI1, outI2=None, qscore=37, fqOffset=33):
    """Grab index sequences for each read in the input fastq and write to separate index fastq files
    Fastq format (per read):
    
    @readName attr1:attr2:attr3:index1[+Index2]
    sequence
    +
    quality score
    
    Output is also fastq with same head line, index1 (or index2) as sequence and a fixed quality score
    """
    quals_index1 = quals_index2 = None # Set to the right length when the first read is processed

    while header := inpFq.stdout.readline():
        name, _, attr = header.strip().partition(" ")
        assert name[0] == "@", "Fastq read name does not start with '@'"
        attrs = attr.split(":")
        assert len(attrs) >= 4, f"Fastq read header does not have 4 fields: {name}"
        index1Seq,_,index2Seq = attrs[3].partition("+")

        if not quals_index1:
            quals_index1 = chr(qscore + fqOffset) * len(index1Seq)
        else:
            assert len(index1Seq) == len(quals_index1), f"Index 1 length is not consistent {name}"
        print(f"{name} {attr}\n{index1Seq}\n+\n{quals_index1}", file=outI1.stdin)

        if outI2:
            assert index2Seq, f"No index2 sequence in fastq: {name}"
            if not quals_index2:
                quals_index2 = chr(qscore + fqOffset) * len(index2Seq)
            else:
                assert len(index2Seq) == len(quals_index2), f"Index 2 length is not consistent: {name}"
            print(f"{name} {attr}\n{index2Seq}\n+\n{quals_index2}", file=outI2.stdin)
        # Ignore remaining lines in the read
        seq = inpFq.stdout.readline()
        plus = inpFq.stdout.readline()
        qual = inpFq.stdout.readline()


def main():
    argparser = argparse.ArgumentParser(description="Extract index sequences from read1 fastq headers and write to separate fastq files")
    argparser.add_argument("read1Fq", help="Input Read1 fastq file", type=pathlib.Path)
    argparser.add_argument("--outDir", help="Output directory", default=".", type=pathlib.Path)
    argparser.add_argument("--no-index2", dest="writeIndex2", help="Don't extract Index2 read", action="store_false")
    args = argparser.parse_args()

    args.outDir.mkdir(parents=True, exist_ok=True)
    assert "_R1" in args.read1Fq.name, "Input should be a read1 fastq (containing '_R1')"
    inpFq = subprocess.Popen(["gzip", "-c", "-d", args.read1Fq], stdout=subprocess.PIPE, text=True)

    indexFq1 = args.outDir / args.read1Fq.name.replace("_R1", "_I1")
    assert not indexFq1.exists(), f"Output file already exists: {indexFq1}"
    outI1 = subprocess.Popen(["gzip", "-c"], stdin=subprocess.PIPE, stdout=open(indexFq1, "w"), text=True)

    if args.writeIndex2:
        indexFq2 = pathlib.Path(args.outDir) / args.read1Fq.name.replace("_R1", "_I2")
        assert not indexFq2.exists(), f"Output file already exists: {indexFq2}"
        outI2 = subprocess.Popen(["gzip", "-c"], stdin=subprocess.PIPE, stdout=open(indexFq2, "w"), text=True)
    else:
        outI2 = None
    writeIndexFqs(inpFq, outI1, outI2)

if __name__ == "__main__":
    main()

