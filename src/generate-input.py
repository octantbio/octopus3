"""
Generates the input fasta "input.fasta" required for full sequence verification.

Does NOT modify original fasta files to be processed.

Appends all fasta files in a specified directory (default is current), 
into a single, sanitized fasta file "input.fasta" which the 
analysis pipeline assumes and hard-codes as input.

Only standard underscores, hyphens, square brackets, and alphanumeric characters 
are allowed in names. Only A, T, C, G, and N are allowed in sequences. 
Invalid characters are replaced with valid characters.
"""

import re
import sys
import argparse
import itertools

from pathlib import Path

# catch broken pipe errors to allow ex) python foo.py ... | head
# see: https://stackoverflow.com/a/30091579
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

#===============================================================================

def fasta_reader(fasta):
    """
    Read in a fasta file lazily and return a generator of the name and sequence
    Input:
    ------
    fasta :: FileType
        opened file

    Yields:
    -------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        read :: str
            Sequence taken from the fasta file

    Requires:
    ---------
    itertools

    Example:
    --------
    itertools.groupby takes a key function and groups all items into a list
    until that key changes. We can key on lines beginning with >, then grab
    every line until the next record in the fasta. This makes our method robust
    to some fasta formats that have forced line breaks at given characters.
    foo = '>ABC>DEF>GHI'
    [(k, list(g)) for k,g in itertools.groupby(foo, lambda x: x == '>')]
    --> [(True, ['>']), (False, ['A', 'B', 'C']), (True, ['>']), ... ]
    Note:
    -----
    Adapted from: https://www.biostars.org/p/710/#1412
    """
    # ditch the boolean (x[0]) and just keep the header/seq grouping
    fa_iter = (x[1] for x in itertools.groupby(fasta, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        name = next(header)[1:].strip()
        # join all sequence lines to one by iterating until the next group.
        read = "".join(s.strip() for s in next(fa_iter))
        yield name, read

#===============================================================================

# silly over optimization to make a fast reverse compliment
# see: https://bioinformatics.stackexchange.com/q/3583
COMP = str.maketrans("ACTGacgt", "TGACtgca")
def rev_comp(seq):
    return seq.translate(COMP)[::-1]

#===============================================================================

if __name__ == '__main__':
    description = 'Generate the input fasta "input.fasta" required for full sequence verification.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("in_dir",
                        nargs='?',
                        default= "./",
                        help="path to directory with all user fastas (default = current directory)")
    parser.add_argument(
        '-o',
        '--out-dir',
        dest='out_dir',
        type=str,
        help='which directory to drop "input.fasta" (default = current directory)',
        default='')
    args = parser.parse_args()

    in_dir = Path(args.in_dir)
    user_files = list(in_dir.glob('*.txt'))
    user_files.extend(list(in_dir.glob('*.fa')))
    user_files.extend(list(in_dir.glob('*.fasta')))
    out_dir = Path(args.out_dir)

    # check if provided directory is valid
    if not in_dir.exists():
        raise FileNotFoundError(f'** provided directory "{str(in_dir)}" does not exist!! **')

    # check if no fastas in provided directory
    if len(user_files) == 0:
        raise Exception("** no user .txt, .fa, or .fasta files found in:" + str(in_dir.absolute()) + "**")

    # checking if input.fasta is already in directory - if so, removing
    if in_dir / "input.fasta" in user_files:
        user_files.remove(in_dir / "input.fasta")

    print("> loaded following files: " + str([str(x).split("/")[-1] for x in user_files]))

    #---------------------------------------------------------------------------

    final_fasta = open(str(out_dir) + "/input.fasta", "w")

    plasmid_names = set()
    plasmid_seqs = {}

    """
    # optional inclusion of specific features of interest to detect for
    # feat is a string of the feature's nucleotide sequence
    # feat_fc is a string of the reverse complement of the feature's nucleotide sequence
    feat = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    feat_fc = "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
    """

    count = 2

    # re-orient the reference plasmids
    for file in user_files:
        fasta = open(file, "r")

        basename = Path(fasta.name).stem
        records = fasta_reader(fasta)

        for name, seq in records:
            seq = seq.upper()
            flat_seq = seq+seq
            # checking for duplicated sequences
            if seq in plasmid_seqs:
                print('> Duplicated Sequences: "' +name+ '" in "' +basename+ '" and "' +plasmid_seqs[seq][0]+ '" in "' +plasmid_seqs[seq][1]+ '" are assigned to the same sequence')
                prompt = input('Would you like to delete"' + name + '" in "' +basename+ '"? (enter y or n): ')
                if prompt in ["yes", "y"]:
                    print("** sequence deleted **")
                    continue
                else:
                    print("** duplicate sequence kept **")
            else:
                plasmid_seqs[seq] = [name, basename]

            # checking for nucleotides in seq other than A, T, C, G, N and replace with N
            (seq, n_sub) = re.subn('[^ATCGN]', 'N', seq)
            if n_sub > 0:
                print('> Invalid Character: "' +name+ '" in "' +basename+ '" contains nucleotide(s) other than A, T, C, G, or N which have been replaced with nucleotide ' + "'N'")

            # checking for spaces in name and replace with underscore '_'
            if " " in name:
                name = name.replace(" ", "_")
                print('> Invalid Character: "' +name+ '" in "' +basename+ '" contains a space which has been replaced with underscore ' + "'_'")

            # checking for parentheses in name and replace with square brackets '[]'
            if "(" in name:
                name = name.replace("(", "[")
                print('> Invalid Character: "' +name+ '" in "' +basename+ '" contains a parenthesis which has been replaced with a square bracket')
            if ")" in name:
                name = name.replace(")", "]")
                print('> Invalid Character: "' +name+ '" in "' +basename+ '" contains a parenthesis which has been replaced with a square bracket')

            # checking for other non-alphanumeric characters and replace with hyphen '-'
            valid = ["-","_","[","]"]
            bad = set()
            for char in name:
                if not char.isalnum():
                    if char not in valid:
                        bad.add(char)
                        name = name.replace(char, "-")
            if len(bad) > 0:
                print('> Invalid Character: "' +name+ '" in "' +basename+ '" contains the following invalid character(s): ' + str(list(bad)) + " which have been replaced with hyphen '-'")

            # checking for duplicated plasmid names
            if name in plasmid_names:
                print('> Duplicate Names: "' +name+ '" in "' +basename+ '" is assigned more than one sequence')
                prompt = input("Would you like to make this name unique? (enter y or n): ")
                if prompt in ["yes", "y"]:
                    name = name + count
                    count += 1
                    plasmid_names.add(name)
                    print("** name changed to: " +name+ " **")
                else:
                    print("** names unmodified **")
            else:
                plasmid_names.add(name)

            """
            # OPTION: check if feat and feat_rc are in sequences and if correctly oriented
            if feat not in flat_seq:
                if feat_rc in flat_seq:  # sequence needs flipped
                    print('> "' +name+ '" in "' +basename+ '" sequence is in wrong direction')
                    prompt = input("Would you like to reverse the sequence? (enter y or n): ")
                    if prompt in ["yes", "y"]:
                        seq = rev_comp(seq)
                        print("** sequnce reversed **")
                    else:
                        print("** sequence unmodified **")
                else:
                    print('> Feature Missing:"' +name+ '" in "' +basename+ '" does not contain feature')
            """

            final_fasta.write(">" + name +"\n")
            final_fasta.write(seq + "\n")

final_fasta.close()
print("*** Complete ***")
print("** input.fasta is located: " + str(out_dir.absolute())+"/input.fasta **")
