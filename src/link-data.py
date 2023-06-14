import re
import sys
import argparse
from pathlib import Path


def clean_name(path):
    """
    Assumes names of form path/to/foo_..._R[1,2]...fastq.*
    """
    clean_path = Path(path)
    head = clean_path.name.split('_')[0]
    read = re.search(r"_(R[1-2])_", clean_path.name).group(1)
    ext = ''.join(clean_path.suffixes)
    return head + "_" + read + ext


if __name__ == '__main__':
    description = 'Create symlinks to all read FastQ\'s and Sample Sheet in a sequencing directory. Just specify where the sequencing run is located relative to the current directory and we\'ll handle the rest!'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('in_dir',
                        type=str,
                        help='path to sequencing run (or stdin if none)')
    parser.add_argument(
        '-o',
        '--out-dir',
        dest='out_dir',
        type=str,
        help='where to drop FastQ and Sample Sheet symlinks (default = current directory)',
        default='')
    args = parser.parse_args()

    # dump links in to a folder with the run id
    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir) / in_dir.name
    out_dir.mkdir(parents=True, exist_ok=True)

    # exit early if out_dir is already organized
    if out_dir.joinpath('organize').exists():
        print(f'** linking skipped since "{str(in_dir)}" is already organized **')
        exit()

    # grab the Sample Sheet
    path_to_samplesheet = in_dir.joinpath('SampleSheet.csv').resolve()
    if path_to_samplesheet.exists():
        try:
            out_dir.joinpath('SampleSheet.csv').symlink_to(path_to_samplesheet)
        except FileExistsError:
            pass
    else:
        raise FileNotFoundError(f'{str(path_to_samplesheet)} does not exist!!')

    # grab the read fastqs
    fastqs = in_dir.glob('Alignment_*/*/Fastq/*.fastq.gz')
    for fq in fastqs:
        out_name = clean_name(fq)
        try:
            out_dir.joinpath(out_name).symlink_to(fq.resolve())
        except FileExistsError:
            pass
