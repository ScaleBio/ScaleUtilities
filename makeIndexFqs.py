#!/bin/env python

import subprocess
import sys
import argparse
import pathlib


def parse_3lvl_header(header):
    """Parse 3lvl format header: @readName attributes(i7+i5)
    Example: @VH02171:15:2227VLHNX:1:1101:19144:1000 1:N:0:GCTCTCGCCT+TCGGATTCGG
    Extracts index sequences from the 4th colon-separated field in attributes
    """
    name, _, attr = header.strip().partition(" ")
    assert name[0] == "@", "Fastq read name does not start with '@'"
    attrs = attr.split(":")
    assert len(attrs) >= 4, f"Fastq read header does not have 4 fields: {name}"
    index1Seq, _, index2Seq = attrs[3].partition("+")
    return name, attr, index1Seq, index2Seq


def parse_qs_header(header):
    """Parse QS format header: @readName:partial-index1 attributes(i7+i5)
    Example: @LH00659:241:22T7WLLT4:1:1101:42065:1140:CTGTCCTAATGGGGTTACCGAAGA 1:N:0:TNCAGACA+GTTCGATA
    Based on OverrideCycles: Y82;I8U24;I8;Y16
    cell-barcode = concatenation of i7 + partial-index1 (from read name)
    i5 = i5 from attributes
    """
    name_part, _, attr = header.strip().partition(" ")
    assert name_part[0] == "@", "Fastq read name does not start with '@'"
    
    # Extract index from read name (after last colon)
    name_parts = name_part.split(":")
    read_name_index = name_parts[-1]  # Last part after colon is the read name index
    
    # Keep the original name with index1 (don't remove it)
    name = name_part
    
    # Extract indices from attributes (4th field)
    attrs = attr.split(":")
    assert len(attrs) >= 4, f"Fastq read header does not have 4 fields: {name}"
    attr_index1, _, attr_index2 = attrs[3].partition("+")
    
    # For QS mode: 
    # index1Seq = cell-barcode
    # index2Seq = i5
    index1Seq = attr_index1 + read_name_index
    index2Seq = attr_index2
    
    return name, attr, index1Seq, index2Seq


def writeIndexFqs_3lvl(inpFq, outI1, outI2=None, qscore=37, fqOffset=33):
    """Process 3lvl format FASTQ files"""
    quals_index1 = quals_index2 = None

    while header := inpFq.stdout.readline():
        name, attr, index1Seq, index2Seq = parse_3lvl_header(header)

        if not quals_index1:
            quals_index1 = chr(qscore + fqOffset) * len(index1Seq)
        else:
            assert len(index1Seq) == len(quals_index1), f"Index 1 length is not consistent {name}"
        
        # Efficient string operation: build complete FASTQ entry
        fastq_entry = f"{name} {attr}\n{index1Seq}\n+\n{quals_index1}\n"
        outI1.stdin.write(fastq_entry)

        if outI2:
            assert index2Seq, f"No index2 sequence in fastq: {name}"
            if not quals_index2:
                quals_index2 = chr(qscore + fqOffset) * len(index2Seq)
            else:
                assert len(index2Seq) == len(quals_index2), f"Index 2 length is not consistent: {name}"
            
            fastq_entry2 = f"{name} {attr}\n{index2Seq}\n+\n{quals_index2}\n"
            outI2.stdin.write(fastq_entry2)
        
        # Ignore remaining lines in the read
        seq = inpFq.stdout.readline()
        plus = inpFq.stdout.readline()
        qual = inpFq.stdout.readline()


def writeIndexFqs_qs(inpFq, outI1, outI2=None, qscore=37, fqOffset=33):
    """Process QS format FASTQ files"""
    quals_index1 = quals_index2 = None

    while header := inpFq.stdout.readline():
        name, attr, index1Seq, index2Seq = parse_qs_header(header)

        # Write index1 (concatenation of attr_index1 + read_name_index)
        if not quals_index1:
            quals_index1 = chr(qscore + fqOffset) * len(index1Seq)
        else:
            assert len(index1Seq) == len(quals_index1), f"Index 1 length is not consistent {name}"
        
        fastq_entry = f"{name} {attr}\n{index1Seq}\n+\n{quals_index1}\n"
        outI1.stdin.write(fastq_entry)

        # Write index2 (attr_index2)
        if outI2:
            assert index2Seq, f"No index2 sequence in fastq: {name}"
            if not quals_index2:
                quals_index2 = chr(qscore + fqOffset) * len(index2Seq)
            else:
                assert len(index2Seq) == len(quals_index2), f"Index 2 length is not consistent: {name}"
            
            fastq_entry2 = f"{name} {attr}\n{index2Seq}\n+\n{quals_index2}\n"
            outI2.stdin.write(fastq_entry2)
        
        # Ignore remaining lines in the read
        seq = inpFq.stdout.readline()
        plus = inpFq.stdout.readline()
        qual = inpFq.stdout.readline()


def writeIndexFqs(inpFq, outI1, outI2=None, qscore=37, fqOffset=33, mode="3lvl"):
    """Main function that routes to appropriate parser based on mode"""
    if mode == "3lvl":
        writeIndexFqs_3lvl(inpFq, outI1, outI2, qscore, fqOffset)
    elif mode == "QS":
        writeIndexFqs_qs(inpFq, outI1, outI2, qscore, fqOffset)
    else:
        raise ValueError(f"Unknown mode: {mode}. Supported modes: 3lvl, QS")


def main():
    argparser = argparse.ArgumentParser(description="Extract index sequences from read1 fastq headers and write to separate fastq files")
    argparser.add_argument("read1Fq", help="Input Read1 fastq file", type=pathlib.Path)
    argparser.add_argument("--outDir", help="Output directory", default=".", type=pathlib.Path)
    argparser.add_argument("--mode", help="FASTQ format mode: 3lvl or QS", choices=["3lvl", "QS"], default="3lvl")
    argparser.add_argument("--no-index2", dest="writeIndex2", help="Don't extract Index2 read", action="store_false")
    args = argparser.parse_args()

    args.outDir.mkdir(parents=True, exist_ok=True)
    
    # Validate input filename based on mode
    if args.mode == "3lvl":
        assert "_R1" in args.read1Fq.name, "Input should be a read1 fastq (containing '_R1')"
    # For QS mode, we don't enforce _R1 requirement as format may vary
    
    inpFq = subprocess.Popen(["gzip", "-c", "-d", args.read1Fq], stdout=subprocess.PIPE, text=True)

    # Create output files based on mode
    if args.mode == "3lvl":
        indexFq1 = args.outDir / args.read1Fq.name.replace("_R1", "_I1")
        if args.writeIndex2:
            indexFq2 = args.outDir / args.read1Fq.name.replace("_R1", "_I2")
    else:  # QS mode
        # Use same logic as 3lvl mode - replace R1 with I1/I2
        indexFq1 = args.outDir / args.read1Fq.name.replace("_R1", "_I1")
        if args.writeIndex2:
            indexFq2 = args.outDir / args.read1Fq.name.replace("_R1", "_I2")

    assert not indexFq1.exists(), f"Output file already exists: {indexFq1}"
    outI1 = subprocess.Popen(["gzip", "-c"], stdin=subprocess.PIPE, stdout=open(indexFq1, "w"), text=True)

    outI2 = None
    if args.writeIndex2:
        assert not indexFq2.exists(), f"Output file already exists: {indexFq2}"
        outI2 = subprocess.Popen(["gzip", "-c"], stdin=subprocess.PIPE, stdout=open(indexFq2, "w"), text=True)

    writeIndexFqs(inpFq, outI1, outI2, mode=args.mode)

if __name__ == "__main__":
    main()
