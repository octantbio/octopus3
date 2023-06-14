![OCTOPUS logo](/img/octopus-logo.png)

[![pytest-integration](https://github.com/octantbio/octopus/actions/workflows/pytest-integration.yml/badge.svg)](https://github.com/octantbio/octopus/actions/workflows/pytest-integration.yml)

(**O**ptimized **C**loning **T**hrough **O**pen-source **P**ipelined **U**nabridged **S**equencing)

## Now introducing: OCTOPUS 3.0
- ONE STEP library prep for easy execution and onboarding
- NEXT-DAY turnaround directly from colonies, faster than Sanger


## Getting Started

- [Installation](/docs/Installation.md)
- [Bench Protocol](/docs/Bench-Protocol.md)
- [Running the Analysis Pipeline](/docs/Running-the-Analysis-Pipeline.md)


## What is OCTOPUS?

OCTOPUS is a light-weight, cost-effective, and robust method for full-plasmid sequence verification using next-generation sequencing.
This respository provides both the bench protocol and complete source code to allow anyone to run this plasmid sequencing pipeline for themselves.

You can read more about the story of OCTOPUS in our [2023 blog post](https://www.octant.bio/blog-posts/octopus-v3), but it's worth summarizing a few features that make OCTOPUS stand out:

- Scalable: one research associate can sequence hundreds—or even thousands—of colonies a week
- Fast: analyzed results in 24 hours from colony/sample submission
- Inexpensive: similar to the cost of a SINGLE Sanger for each sample
- Robust: our RCA-based protocol performs exceedingly well and automated QC can detect common colony-picking errors
- Easy: OCTOPUS can be performed at scale with common lab equipment (no upfront capital investment in esoteric instruments)
- Straighforward: OCTOPUS works directly off of picked colonies and crude lysates (no automated DNA purification required)

Did we mention scalable? OCTOPUS has allowed us to sequence thousands of colonies a month for now over four years, and usage is still strong!

![OCTOPUS runs over time 3.0](/img/cumulative-wells-sequenced.png)

### How OCTOPUS Works

![OCTOPUS 3.0 timeline vs v2.0](/img/run-timeline-compare.png)

OCTOPUS 3.0 consists of the following steps:
1. Colony picking and glycerol stock setup
2. Fast RCA
3. seqWell [ExpressPlex](https://seqwell.com/expressplex-library-prep-kit/) 90-minute library prep
4. Illumina sequencing
5. Turn-key data analysis
6. Pick perfect sequences _next day_ from colonies

Interested in learning more?
[Check out our documentation for more information](/docs/).

## Contributing

We welcome any feedback to help make OCTOPUS better!
If you encounter problems while running this pipeline or have suggestions for further improvements, please open an issue using our [issue tracker](https://github.com/octantbio/octopus3/issues).

Pull requests are also welcome.

## License

This project (including the bench protocol and source code) is licensed under the Apache 2.0 License - see the [LICENSE](/LICENSE) file for details. Additional licensing information:

- [BBTools](docker/bbtools-license) - custom
- [mlr](docker/mlr-license) - custom
- [starcode](docker/starcode-license) - GPL-3.0
- [SPAdes](docker/spades-license) - GPL-2.0
- [minimap2](https://github.com/lh3/minimap2/blob/master/LICENSE.txt) - MIT
- [samtools](https://github.com/samtools/samtools/blob/develop/LICENSE) - MIT
- [bcftools](https://github.com/samtools/bcftools/blob/develop/LICENSE) - MIT
- [htslib](https://github.com/samtools/htslib/blob/develop/LICENSE) - MIT
- [freebayes](https://github.com/ekg/freebayes/blob/master/LICENSE) - MIT


[OCTOPUS logo](/img/octopus-logo.png) copyright Octant Bio. All rights reserved.

