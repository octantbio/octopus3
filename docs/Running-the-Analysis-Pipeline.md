**Running the Analysis Pipeline**
=================================

We provide an analysis pipeline that can get you kickstarted in finding the perfect clones!

## **Method 1: run-octopus-analysis.sh**

We can walk through this method with an example:

Inside the provided `test/` directory happens to contain the two necessary components for a run with run-id `example_octopus2_miseq_run`:

1. The sequencing data located in `/test/miseq-outputs`
2. Some possible reference sequences in `/test/fastas`

To run analysis for this run, follow the following steps:

1. `cd` to the root of the `octopus3` directory if not yet already.
    ```
    cd /path/to/octopus3
    ```
2. Assuming the docker image name happens to be `octopus3:release`, run the script with the full absolute directories to the sequencing data and reference sequences.
    ```
    ./run-octopus-analysis.sh -a $(pwd)/test/fastas/example_octopus2_miseq_run -q $(pwd)/test/miseq-outputs/example_octopus2_miseq_run -d octopus3:release
    ```
    - Run `./run-octopus-analysis.sh -h` to see a help screen for more details
    - Not providing the fasta location will run _de novo_ assembly only.
3. Wait for analysis and then done! Obtain the output `aggregated-stats.tsv` in `/pipeline/example_octopus2_miseq_run` directory.

## **Method 2: Manually**

### **1. Link Sequencing Data**

As OCTOPUS uses `make` to orchestrate everything, there are some conventions your data must adhere to. First, we need to link the output folder from your sequencing run to the `pipeline/` directory. This can be done with `/src/`[link-data.py](/src/link-data.py), which when provided the root directory of the sequencing data, will create symlinks to `SampleSheet.csv` and `*.fastq.gz` files in `/pipeline/<run-sequencing-folder>`. Let's use the `test/` provided examples from [Method 1](#method-1-run-octopus-analysissh) again.

```
cd /path/to/octopus3

python3 src/link-data.py test/miseq-outputs/example_octopus2_miseq_run -o pipeline/
```

In general, we typically use the default folder name produced by the sequencer as an identifier. Second, many steps in the OCTOPUS pipeline will process the file name of the fastq's. So to avoid issues, make sure all unique information is contained before the first underscore in your SampleSheet (most Illumina sequencers will automatically convert any _'s in the `Sample_Name` column of the SampleSheet to -'s anyways). Importantly, the pipeline will trim out anything between the first underscore and the read specifier (e.g. `my-reads_foo_bar_baz_R1.fastq.gz -> my-reads_R1.fastq.gz`) to ensure everything behaves properly.

### **2. Generate Reference Library**

Next, locate the directory containing the reference sequences in fasta format of the plasmids you are trying to sequence. These fastas need to be sanitized for invalid characters, duplicates, and misassignments as well as be combined into a single fasta called `input.fasta` placed in `/pipeline/<run-sequencing-folder>` where the symlinks from earlier also reside. This can all be done with `/src/`[generate-input.py](/src/generate-input.py).

```
python3 src/generate-input.py test/fastas/example_octopus2_miseq_run -o pipeline/example_octopus2_miseq_run
```

Alternatively, you can manually place `input.fasta` in `/pipeline/<run-sequencing-folder>`.

### **_De Novo_ Assembly**

If you do not know or have your input, run `make denovo` instead of `make all` in the later step to take the pipeline through the _de novo_ assembly step only.

### **3. Running the Pipeline**

After getting the data in place, make note again of the path to the root of the sequencing data (you had this when [linking data in step 1](#1-link-sequencing-data)) but this time its full absolute path.

From there, we can drop into the docker image [that was built before all these steps](/docs/Installation.md). Let's assume it has the name `octopus3:release`.

```
cd /path/to/octopus3

docker run --rm -it -v $(pwd)/test/miseq-outputs/example_octopus2_miseq_run:$(pwd)/test/miseq-outputs/example_octopus2_miseq_run -v $(pwd):/data/octopus3 octopus3:release bash
```

This links your octopus folder (`/path/to/octopus3`) to the docker image (`/data/octopus`) as well as ensures the path to the sequencing data is exactly preserved so that the [symlinks generated in step 1](#1-link-sequencing-data) are not broken. Note that docker requires you to specify the *absolute* path to the directory (`$(pwd)` is a handy shortcut to do that for you). Also, be aware that `--rm` makes the image ephemeral so anything written outside of the octopus directory will be lost if you log out of the container. From the docker image, we can

```
cd octopus3

make all
```

to run the pipeline on every sequencing run in the `pipeline/` directory, creating the output for the `test/` run as `/pipeline/example_octopus2_miseq_run/aggregated-stats.tsv`.

## **Output: aggregated-stats.tsv**

As the name suggests, results pertinent to an OCTOPUS run are aggregated into a `tsv` file for your analysis. The columns are:

- `Run`: Illumina run ID
- `Plate`: plate ID
- `Well`: well address
- `Plate_Well`: unique plate_well identifier for each sample
- `DeNovo_Ref`: well identity based on aligning _de novo_ assembly to reference library
- `CIGAR`: CIGAR string from aligning the _de novo_ assembly to `DeNovo_Ref`
- `LT_10`: percentage of input reference with < 10x coverage (ideally close to 0)
- `LT_3`: percentage of input reference sequence with < 3x coverage (if not 0 inspect read pileup)
- `BC_Contam`: are there multiple barcodes in this well (TRUE/FALSE)? ([more details](/docs/Pipeline-Details.md#barcode-filter))
- `n_vars`: number of variants detected by FreeBayes (note barcodes count as variants)
- `n_barcodes`: number of barcodes detected
- `expected_bcs`: expected number of barcodes based on the reference (in a perfect plasmid `n_vars == n_barcodes == expected_bcs`)
- `bc_1`: sequence of barcode 1 pulled from the variant caller (may be reverse complement)
- `pos_1`: position of barcode 1 in _de novo_ assembly
- `bc_N`: sequence of barcode N pulled from the variant caller (may be reverse complement; NA if missing)
- `pos_N`: position of barcode N in _de novo_ assembly (NA if missing)
- `Contaminants`: number of reads from "contaminants" ([more details](/docs/Pipeline-Details.md#alternative-contaminants))
- `Leftover`: number of reads in well leftover after filtering out "contaminants"
- `Consensus`: the reference sequence but with detected variants substituted in
- `Contig`: the _de novo_ assembly. Note the first and last N bases (often 55 or 125) are repeated
- `Ref_Seq`: sequence that _de novo_ assembly aligns to, aka the reference sequence

## **Picking perfects**

One way you can analyze the results is by pasting the `aggregated-stats.tsv` into a spreadsheet

1. If applicable, filter out any "TRUE" or "NA" values under `BC_Contam`
2. If applicable, flag or filter out any duplicate barcodes
3. Filter out any unexpected variants. The pipeline will automatically detect any strings of N's in the `input.fasta` and report the number of `expected_bcs` for that reference. Perfect clones should have `expected_bcs == n_barcodes == n_vars`
4. Ensure that there is adequate coverage by checking `LT_10` and `LT_3`. We recommend only picking wells with `LT_3 == 0`. You can be more conservative by using `LT_10` to specify your cutoff. (For a 10kb plasmid an `LT_3` of 0.001 means that 10 bp of the plasmid did not have a coverage of at least three).
5. If there happens to be a clone that does not have sufficient coverage (by `LT_10` or `LT_3`) but is absolutely required, use `samtools tview` to manually inspect the read pileup in critical areas of your plasmid:
    1. In a new terminal, `cd` into your octopus directory
    2. Open up a new docker instance - `docker run --rm -it -v "$(pwd)":/root/octopus octant/octopus3 /bin/bash`
    3. Navigate into the folder that contains the analyzed data from that run - `cd octopus/pipeline/your_run_id`
    4. View the pileup - `samtools tview data/your_plate/your_well.map.bam lib/your_ref.fasta`
        - Note `your_ref` will be `DeNovo_Ref` in `aggregated-stats.tsv`

