#!/usr/bin/env python3
# Task: Replace all instances of 'a,c,g,t' in a FASTA file with 'A,C,G,T'

import sys

def upcase_sequence():
    """Takes standard input, transforms sequence characters to uppercase.

    Use with pipes on UNIX commandline as a filter to transform lowercase
    sequence to uppercase sequence.
    """
    for line in sys.stdin:
        if not line.startswith(">"):
            line = line.upper()
        sys.stdout.write(line)

if __name__ == '__main__':
    upcase_sequence()
