**OCTOPUS 3.0 Bench Protocol**
==============================

# Table of Contents

- [Reagents & Consumables](#reagents-and-consumables)
- [Equipment](#equipment)
- [Day 1 - Input, Library, and Sequencing Prep](#day-1---input-library-and-sequencing-prep)
    - [1. Input Generation: Colony Picking and Fast RCA](#1-input-generation-colony-picking-and-fast-rca)
        - [1.1 Colony preparation](#11-colony-preparation)
        - [1.2 RCA reaction](#12-rca-reaction)
    - [2. NGS Library Prep: seqWell ExpressPlex Kit](#2-ngs-library-prep-seqwell-expressplex-kit)
        - [2.1 ExpressPlex setup](#21-expressplex-setup)
        - [2.2 Pooling](#22-pooling)
    - [3. Library Quantification and Sequencing Setup](#3-library-quantification-and-sequencing-setup)
- [Day 2 - Run Wrapup](#day-2---run-wrapup)
    - [4. Glycerol Stock Culture Plates](#4-glycerol-stock-culture-plates)
    - [5. Run Analysis Pipeline](#5-run-analysis-pipeline)

# Reagents and Consumables:
The following are just where we happen to source our materials. Some specific plate types are assumed by the OT-2 protocol we provide.

- [seqWell ExpressPlex Library Prep Kit](https://seqwell.com/expressplex-library-prep-kit/) from seqWell
- [Thermo Fisher EquiPhi29 DNA Polymerase](https://assets.thermofisher.com/TFS-Assets%2FLSG%2Fmanuals%2FMAN0017945_EquiPhi29_DNA_Polymerase_PI.pdf) kit from Thermo Fisher
- [Illumina MiSeq Reagent Micro Kit v2 (300-cycles)](https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/miseq-reagent-kit-v2.html) (MS-103-1002) sufficient for up to 8x96 samples
- [Biosearch 25 mM dNTP set](https://shop.biosearchtech.com/pcr-instruments%2c-reagents-and-consumables/pcr-and-qpcr-reagents/pcr-and-qpcr-enzymes-and-reagents/premixed-dntp-solutions/p/D59104)
- [MCLAB 500 µM Exo Resistant Random Hexamers](https://www.mclab.com/Exo-Resistant-Random-Hexamer.html)
    - We found that having a small amount of random hexamers normalizes against amplification bias of regions of plasmid targeted by the custom RCA primers and increases yield.
- [IDT 100µM custom RCA primers](/docs/Focused-RCA.md) (OPTIONAL): Your own custom "RCA Primers" mixed equimolarly each at 100uM, designed to bind to a common sequence on your plasmids
    - For example, we designed 8 primers that target the E1 origin of replication that is in all our plasmids.
    - We ordered these from [IDT](https://www.idtdna.com/pages/support/faqs/how-do-i-order-a-phosphorothioate-modified-oligo-) with the last three 3' bases phosphorothiorated to make them resistant to phi29 3'->5' exonuclease activity.
- [Flat bottom 96-well plates](https://www.sigmaaldrich.com/catalog/product/aldrich/br781602)
- [Bio-Rad Hard-Shell 96-well PCR plates](https://www.bio-rad.com/en-us/sku/HSP9601-hard-shell-96-well-pcr-plates-low-profile-thin-wall-skirted-white-clear)
- [Applied Biosystems MicroAmp EnduraPlate 384-well PCR plates](https://www.thermofisher.com/order/catalog/product/A36931)
- [Polarseal foil plate seal](https://www.thomassci.com/Molecular-Diagnostics/Amplification/PCR/Plate-Sealing/_/PolarSeal-Aluminum-Microplate-Seals) for thermocycling
- [Denovix High Sensitivity dsDNA Assay kit](https://www.denovix.com/purchase-fluorescence-assays/#1600940711334-569803f2-3d7a)
- [NEST 2 mL 96-well deep well plate, V-bottom](https://shop.opentrons.com/nest-2-ml-96-well-deep-well-plate-v-bottom/) (if pooling by OT-2)
- PCR strip tubes (optional)
---
- ddH2O
- Ethanol
- [IDTE (10 mM Tris, 0.1 mM EDTA, pH 8.0)](https://www.idtdna.com/pages/products/reagents-and-kits/buffers-and-solutions)
- 50% [Glycerol](https://www.fishersci.com/shop/products/glycerol-molecular-biology-fisher-bioreagents-2/BP2291)
- [2xYT](https://us.vwr.com/store/product/7437420/vwr-life-science-2xyt-medium-broth)
- 1 M [Sodium Hydroxide](https://www.fishersci.com/shop/products/sodium-hydroxide-pellets-certified-acs-fisher-chemical-7/S318100)

# Equipment:

- [96-well thermocycler](https://www.bio-rad.com/en-us/product/c1000-touch-thermal-cycler?ID=LGTW9415)
- [384-well thermocycler](https://www.thermofisher.com/order/catalog/product/4388444)
- 37˚C shaking incubator
- Nanodrop, Qubit, BioAnalyzer or Tapestation
- [Illumina MiSeq](https://www.illumina.com/systems/sequencing-platforms/miseq.html)
- Plate centrifuge
- Vortex mixer
- Multichannel 10
- Multichannel 200
- [Rainin Liquidator 20](https://www.shoprainin.com/Products/Pipettes-and-Tips/Pipettes/High-throughput-Pipetting/Liquidator%E2%84%A2-96/Liquidator-96-2-20-%C2%B5L-LIQ-96-20/p/17014207) (optional)
- [Rainin Liquidator 200](https://www.shoprainin.com/Products/Pipettes-and-Tips/Pipettes/High-throughput-Pipetting/Liquidator%E2%84%A2-96/Liquidator96%2C-5-200-%C2%B5L-LIQ-96-200/p/17010335) (optional)
- [Integra Viaflo 384](https://shop.integra-biosciences.com/us/s/product/detail/01tD0000005XbLCIA0) (optional)
- [Opentrons OT-2](https://shop.opentrons.com/ot-2-robot/) (optional)
- Magnetic stand for 1.5 ml tubes (optional)

# Day 1 - Input, Library, and Sequencing Prep
## 1. Input Generation: Colony Picking and Fast RCA

### **1.1 Colony preparation**
1. Prepare "Sample Plate": Fill each well of a 96-well PCR plate with 25 µl of water. This will be the input to OCTOPUS.
2. Prepare "Culture Plate": Fill each well of a 96-well flat-bottom plate with 110 µl of culture media + antibiotic (we use 2xYT). This will be for glycerol stock storage.
3. Pick colonies into Culture Plate, mix well, then transfer 5 µl of culture into the Sample Plate.
4. Grow Culture Plate overnight shaking at 30-37˚C with lid taped on to prevent evaporation.
5. Seal Sample Plate and lyse samples by heating at 95˚C for 3 minutes.

---
Notes:
- The extra 10 µl culture media is to roughly compensate for evaporation.
- We cover each Culture Plate with its lid and tape tightly down the middle of all four edges of the plate. We've tried breathable seals at one point but found noticeable evaporation, but YMMV.
- Plates are set aside at RT in a safe area with low risk of bumping before being put into a large shaking incubator to saturate overnight.
---

### **1.2 RCA reaction**
1. Prepare 200 µl 2X RCA reaction for each 96-well plate:

<table>
<tr>
<th> Normal RCA </th>
<th> "Focused" RCA </th>
</tr>
<tr>
<td>

- 88 µl water
- 40 µl EquiPhi29 10X Reaction Buffer
- 40 µl 10 mM dNTP mix
- 4 µl 100 mM DTT
- 8 µl 500 µM Exo-Resistant Random Hexamers
- 20 µl EquiPhi29 DNA Polymerase (0.1 µg/µl)

</td>
<td>

- 93.4 µl water
- 40 µl EquiPhi29 10X Reaction Buffer
- 40 µl 10 mM dNTP mix
- 4 µl 100 mM DTT
- 1.6 µl 100 µM custom RCA Primers
- 1 µl 500 µM Exo-Resistant Random Hexamers
- 20 µl EquiPhi29 DNA Polymerase (0.1 µg/µl)

</td>
</tr>
</table>

2. Add 2 µl of the 2X RCA reaction mixture to each well of a 384-well plate. Centrifuge briefly.
3. Add 2 µl of the lysed sample to each well of the "RCA Plate," now a total 4 µl reaction volume. Centrifuge briefly.
4. Incubate RCA Plate on thermocycler at 42˚C for 3 hours, then heat-inactivate at 65˚C for 10 minutes.

---
Notes:
- We use the Rainin Liquidator for 96-well transfers (e.g. Sample Plate -> RCA Plate).
- RCA Plate can be pre-made and stocked at -20˚C, stable for at least a month but likely a lot longer. Make sure to spin down pre-stocked plates well after thawing to make sure mixture is at the bottom of every well.
- We've only vetted the Thermo Fisher EquiPhi29 RCA kit for fast RCA so far. It is miniaturized 5X (20 µl to 4 µl) and modified to contain, in addition to the random hexamers, "focused" primers specific to the ColE1 origin of replication that allow for enhanced amplification of plasmid DNA as opposed to gDNA or other contaminants. See the [Focused RCA](/docs/Focused-RCA.md) page for more details.
- If not running a multiple of 384 samples, fill the 384-well plate in the following pattern: Quadrant A1, A2, B1. This format allows for plate-wide transfers and pooling by OT-2.
    <details>
    <summary><strong>96-to-384 plate pattern diagram</strong></summary>

    ![96-to-384 plate pattern](/img/quadrants-384-well-plate.png)
    [diagram source](https://hwpi.harvard.edu/files/iccb/files/384-to-96_library_plate_format.png)
    </details>

    We have not tested these miniaturized reactions in the 96-well plate, only the 384.
---

## 2. NGS Library Prep: seqWell ExpressPlex Kit

### **2.1 ExpressPlex setup**
1. Gently add 2 µl of Ready Reaction Plate (red) into a 384-well PCR "Reaction Plate." Be gentle, as reagent easily bubbles.
2. Gently add 1 µl of desired Indexing Plate (white) into corresponding wells of Reaction Plate.
3. Make a 1:16 dilution of RCA Plate, recommended in water. [Warning: RCA Plate wells may be very goopy.]
4. Add 1 µl from the RCA 1:16 dilution plate into Reaction Plate, now a total 4 µl reaction volume. Centrifuge briefly.
5. Run Reaction Plate on thermocycler with the ExpressPlex protocol:
    | Step # | Temp (°C) | Time (mm:ss) | # of Cycles |
    | ------ | --------- | ------------ | ----------- |
    | 1      | 55        | 5:00         | 1           |
    | 2      | 75        | 5:00         | 1           |
    | 3      | 79        | 5:00         | 1           |
    | 4      | 83        | 5:00         | 1           |
    | 5      | 98        | 3:00         | 1           |
    | 6      | 98        | 0:15         | 12          |
    |        | 64        | 0:30         |             |
    |        | 72        | 1:00         |             |
    | 7      | 72        | 5:00         | 1           |
    | 8      | 4         | hold         | -           |

---
Notes:
- Protocol is adapted from the [ExpressPlex protocol](https://5527966.app.netsuite.com/core/media/media.nl?id=242291&c=5527966&h=OBqIoVnIJtTd8LUnt_rQUEUHw-385wvvAxRUzpb0p3QAOE7-&_xt=.pdf).
- ExpressPlex reaction is miniaturized 4X (16 µl to 4 µl).
    - We have not tested these miniaturized reactions in the 96-well plate, only the 384.
- We use the Rainin Liquidator for 96-well transfers (e.g. Ready Reaction Plate, Indexing Plate -> Reaction Plate).
- We use the Integra Viaflo for 384-well transfers (e.g. RCA Plate -> dilution plate -> Reaction Plate).
- The 1:16 RCA dilution typically results in a final DNA concentration of ~5-15 ng/µl for us. Any input <4 ng/µl we've found to be less reliable.
- To dilute RCA input, we thoroughly mix 1 µl from RCA Plate into a separate 384-well PCR plate containing 15 µl IDTE. IDTE dilution happens to be more convenient for our workflow and hasn't been a problem for us but YMMV, as seqWell doesn't officially recommend diluting input in EDTA-containing solutions.
---

### **2.2 Pooling**
**Method 1: Automated by OT-2**

The [Magnetic Module](https://shop.opentrons.com/magnetic-module-gen2/) is required.

See the [OT-2 Setup](/docs/OT-2-Setup.md) page for details.

**Method 2: By hand**
1. Bring MAGwise beads (included in ExpressPlex kit) to room temperature if not already. Fully resuspend by vortexing.
2. Use multichannel to consolidate all contents of each indexed 96 wells of Reaction Plate into their own labeled PCR strip tube (i.e. 1 strip tube per plate index). This will be bubbly.
3. Further consolidate each strip tube and bubbles into a labeled 1.5 ml LoBind tube. Centrifuge thoroughly to collapse excess bubbles.
4. Ensure MAGwise beads are fully resuspended. Add 0.75 volumetric equivalent of beads to each tube (e.g. 225 µl beads into 300 µl pool). Pipette mix thoroughly and incubate on bench for 5 minutes.
5. Place tubes on magnetic stand for 5 minutes to pellet.
6. Keeping on stand, discard supernatant without disturbing pellet.
7. Keeping on stand, add 80% ethanol over pellet enough to cover it. After ≥30 seconds, discard supernatant without disturbing pellet.
8. Repeat the previous ethanol wash step for a total of 2 washes.
9. Quick spin tubes and place back on stand. Use small pipette tip to remove residual ethanol. Immediately proceed to next step.
10. Take tubes off stand, thoroughly resuspend pellet with pipette in 30 µl IDTE, and incubate for 5 minutes to elute.
11. Place tubes on stand for 2 minutes and transfer supernatant to new labeled tubes. These are the final libraries.

---
Notes:
- If this is your first time performing OCTOPUS, we recommend running gel samples after step 3 and step 11. You should observe a smear concentrated at around 400-600 bp.
---

## 3. Library Quantification and Sequencing Setup

To prepare the libraries for loading the MiSeq, we follow Protocol A of [Illumina's Denature and Dilute Libraries Guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/miseq-denature-dilute-libraries-guide-15039740-10.pdf) for the v2 workflow and 4 nM dilution. We aim for a cluster density of ~800 K/mm2 and spike in 30 µl denatured 12.5 pM PhiX to the final 600 µl diluted library.

The sequencing run takes ~18.5 hours with the [Illumina MiSeq Reagent Micro Kit v2 (300-cycles)](https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/miseq-reagent-kit-v2.html). We found that we can reliably load up to 8x96 samples per run with the Micro flow cell. We max out the [329 cycles available with the 300-cycle kit](https://knowledge.illumina.com/instrumentation/general/instrumentation-general-reference_material-list/000007002):
- 151x2 cycles to 151 bp paired-end reads
- 10 cycles to i7 index read
- 10 cycles to i5 index read
- 7 extra cycles needed for dual indexing with this kit

We found that variability in quantification can be reduced by quantifying a previously run OCTOPUS library and using that as a baseline to estimate the concentration of the current libraries. In our hands, a fluorescence based assay like the [Denovix High Sensitivity dsDNA Assay kit](https://www.denovix.com/purchase-fluorescence-assays/#1600940711334-569803f2-3d7a) gives sufficiently accurate quantification alongside a previously run library. If no previously run libraries are available, quantification by qPCR has been the most accurate method. Every quantification method is slightly different and each comes with its own pros and cons. See [Illumina's guide on best practices for library quantification](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000003750).

The Sample Sheet contains an entry for every combination of 96-well and plate-wide index. For each run, we 
1. Make a copy of an older Sample Sheet
2. Update the ```[Header]``` section fields (e.g. Experiment Name, Date, etc.)
3. Find and replace the older plate names in the ```[Data]``` section with the current ones
4. Find and replace the ```I5_Index_ID``` and ```index2``` fields of the ```[Data]``` section with the ones used for the current run
    -  See the Resources section of the [ExpressPlex product page](https://seqwell.com/products/expressplex-library-prep-kit/) for the most up-to-date master list of indexes.
5. The ```I7_Index_ID``` and ```index``` fields of the ```[Data]``` section should be constant for every 96 set of samples

Feel free to remove any fields you find redundant/unnecessary, though keep at least the ```Sample_Plate``` field if using the analysis pipeline we provide since it relies on extracting info from that field to organize the FastQs by plate.

See the dropdown below for a template Sample Sheet that contains 2x96 samples:

<details>
<summary><strong>Example SampleSheet.csv with 2x96 samples</strong></summary>

```
[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Investigator Name,Octonaut,,,,,,,,
Project Name,Misc,,,,,,,,
Experiment Name,20230613_OCTOPUS_plate001-002,,,,,,,,
Date,6/13/2023,,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,Nextera,,,,,,,,
Description,,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
151,,,,,,,,,
151,,,,,,,,,
[Settings],,,,,,,,,
Adapter,CTGTCTCTTATACACATCT,,,,,,,,
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
OCTOPUS_plate001_A01,plate001_A01,plate001,A01,i7-1437,GTCAAGTCCA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B01,plate001_B01,plate001,B01,i7-1610,TATCTCTTCC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C01,plate001_C01,plate001,C01,i7-0790,CCGCGAAGAA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D01,plate001_D01,plate001,D01,i7-0822,CCTACTCGGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E01,plate001_E01,plate001,E01,i7-1872,TTCGTATCAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F01,plate001_F01,plate001,F01,i7-1666,TCCTCCATCC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G01,plate001_G01,plate001,G01,i7-0474,ATAACATCGC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H01,plate001_H01,plate001,H01,i7-1201,GATATGCGTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A02,plate001_A02,plate001,A02,i7-0598,CAACTAACTC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B02,plate001_B02,plate001,B02,i7-1419,GTACTGGATT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C02,plate001_C02,plate001,C02,i7-0708,CATCGGAGGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D02,plate001_D02,plate001,D02,i7-0574,ATTCTGATGG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E02,plate001_E02,plate001,E02,i7-1608,TATCGTTACC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F02,plate001_F02,plate001,F02,i7-1789,TGGTCTGTTA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G02,plate001_G02,plate001,G02,i7-1758,TGCCAACATG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H02,plate001_H02,plate001,H02,i7-1643,TCATTACACG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A03,plate001_A03,plate001,A03,i7-0475,ATAACCTGAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B03,plate001_B03,plate001,B03,i7-1768,TGCGGTTCCA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C03,plate001_C03,plate001,C03,i7-0127,AATACTTGCC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D03,plate001_D03,plate001,D03,i7-0139,AATCGCGGAA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E03,plate001_E03,plate001,E03,i7-0661,CAGAACGCGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F03,plate001_F03,plate001,F03,i7-1722,TGAACCAAGG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G03,plate001_G03,plate001,G03,i7-1819,TGTTAGTCAG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H03,plate001_H03,plate001,H03,i7-0683,CAGTAGGTAA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A04,plate001_A04,plate001,A04,i7-0680,CAGGTACTTC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B04,plate001_B04,plate001,B04,i7-0678,CAGGAATATG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C04,plate001_C04,plate001,C04,i7-0043,AACGCACAAT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D04,plate001_D04,plate001,D04,i7-0576,ATTGAGAAGG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E04,plate001_E04,plate001,E04,i7-0578,ATTGCACCTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F04,plate001_F04,plate001,F04,i7-0116,AAGTGGATAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G04,plate001_G04,plate001,G04,i7-1192,GATAAGATGC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H04,plate001_H04,plate001,H04,i7-1555,TACCTCGACA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A05,plate001_A05,plate001,A05,i7-0028,AACCGAGCCA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B05,plate001_B05,plate001,B05,i7-1782,TGGATTCAAG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C05,plate001_C05,plate001,C05,i7-0586,CAACAGATAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D05,plate001_D05,plate001,D05,i7-0102,AAGGTAACTC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E05,plate001_E05,plate001,E05,i7-1728,TGACAATACG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F05,plate001_F05,plate001,F05,i7-1221,GATTGTGCAT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G05,plate001_G05,plate001,G05,i7-1231,GCACACCATT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H05,plate001_H05,plate001,H05,i7-1695,TCGTTATTCC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A06,plate001_A06,plate001,A06,i7-0597,CAACGTCATT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B06,plate001_B06,plate001,B06,i7-0068,AAGAACGATG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C06,plate001_C06,plate001,C06,i7-1059,CTGCAATTAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D06,plate001_D06,plate001,D06,i7-0907,CGCCGATGAT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E06,plate001_E06,plate001,E06,i7-0928,CGGACAAGAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F06,plate001_F06,plate001,F06,i7-0952,CGTACTCCTC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G06,plate001_G06,plate001,G06,i7-1757,TGCATGAGTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H06,plate001_H06,plate001,H06,i7-0141,AATCTGGAGC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A07,plate001_A07,plate001,A07,i7-0581,ATTGGTCAGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B07,plate001_B07,plate001,B07,i7-0075,AAGACCTGTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C07,plate001_C07,plate001,C07,i7-0917,CGCTAATGAA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D07,plate001_D07,plate001,D07,i7-0668,CAGAGTGCAT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E07,plate001_E07,plate001,E07,i7-1364,GGTAAGCTGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F07,plate001_F07,plate001,F07,i7-1432,GTATTCAGTG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G07,plate001_G07,plate001,G07,i7-1286,GCGCTACGTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H07,plate001_H07,plate001,H07,i7-1900,TTGTCATAGC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A08,plate001_A08,plate001,A08,i7-0558,ATGTTGCGGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B08,plate001_B08,plate001,B08,i7-0147,AATGCTAACC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C08,plate001_C08,plate001,C08,i7-1404,GTAACACGTA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D08,plate001_D08,plate001,D08,i7-1199,GATATACGGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E08,plate001_E08,plate001,E08,i7-0513,ATCCGAGAGG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F08,plate001_F08,plate001,F08,i7-0623,CAATTCACAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G08,plate001_G08,plate001,G08,i7-1899,TTGTCAGTTC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H08,plate001_H08,plate001,H08,i7-0485,ATACGCCATT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A09,plate001_A09,plate001,A09,i7-0636,CACCAATAAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B09,plate001_B09,plate001,B09,i7-1296,GCGTCCACAA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C09,plate001_C09,plate001,C09,i7-0555,ATGTGCGCTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D09,plate001_D09,plate001,D09,i7-1609,TATCTAGTGC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E09,plate001_E09,plate001,E09,i7-1268,GCCTACAATG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F09,plate001_F09,plate001,F09,i7-0096,AAGCTCAGTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G09,plate001_G09,plate001,G09,i7-1562,TACGCTTAGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H09,plate001_H09,plate001,H09,i7-0510,ATCCACTAGG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A10,plate001_A10,plate001,A10,i7-1802,TGTCCGTCTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B10,plate001_B10,plate001,B10,i7-0716,CATGAGTAAC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C10,plate001_C10,plate001,C10,i7-1063,CTGCGCGAAT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D10,plate001_D10,plate001,D10,i7-1074,CTGTAGTATG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E10,plate001_E10,plate001,E10,i7-1173,GAGCCGTACA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F10,plate001_F10,plate001,F10,i7-0727,CATTCTTAGG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G10,plate001_G10,plate001,G10,i7-1730,TGACCGACAA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H10,plate001_H10,plate001,H10,i7-1816,TGTGTAACCG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A11,plate001_A11,plate001,A11,i7-0862,CGAAGGACTG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B11,plate001_B11,plate001,B11,i7-1696,TCTACCGTCA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C11,plate001_C11,plate001,C11,i7-1365,GGTAATATCG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D11,plate001_D11,plate001,D11,i7-1813,TGTGCGAGTT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E11,plate001_E11,plate001,E11,i7-1848,TTATCGCTGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F11,plate001_F11,plate001,F11,i7-1400,GGTTGAGTTC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G11,plate001_G11,plate001,G11,i7-0580,ATTGGACGCC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H11,plate001_H11,plate001,H11,i7-1616,TATGTGTGTG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_A12,plate001_A12,plate001,A12,i7-1598,TAGTTATCGC,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_B12,plate001_B12,plate001,B12,i7-1755,TGCAGGTGAT,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_C12,plate001_C12,plate001,C12,i7-1809,TGTGAAGCTA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_D12,plate001_D12,plate001,D12,i7-0599,CAACTCCTGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_E12,plate001_E12,plate001,E12,i7-1795,TGTAGCAACG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_F12,plate001_F12,plate001,F12,i7-0937,CGGTAACGCA,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_G12,plate001_G12,plate001,G12,i7-1423,GTAGCAGCAG,i5-1714,GTAACACAGA,,
OCTOPUS_plate001_H12,plate001_H12,plate001,H12,i7-0989,CTACAGCCGA,i5-1714,GTAACACAGA,,
OCTOPUS_plate002_A01,plate002_A01,plate002,A01,i7-1437,GTCAAGTCCA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B01,plate002_B01,plate002,B01,i7-1610,TATCTCTTCC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C01,plate002_C01,plate002,C01,i7-0790,CCGCGAAGAA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D01,plate002_D01,plate002,D01,i7-0822,CCTACTCGGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E01,plate002_E01,plate002,E01,i7-1872,TTCGTATCAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F01,plate002_F01,plate002,F01,i7-1666,TCCTCCATCC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G01,plate002_G01,plate002,G01,i7-0474,ATAACATCGC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H01,plate002_H01,plate002,H01,i7-1201,GATATGCGTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A02,plate002_A02,plate002,A02,i7-0598,CAACTAACTC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B02,plate002_B02,plate002,B02,i7-1419,GTACTGGATT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C02,plate002_C02,plate002,C02,i7-0708,CATCGGAGGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D02,plate002_D02,plate002,D02,i7-0574,ATTCTGATGG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E02,plate002_E02,plate002,E02,i7-1608,TATCGTTACC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F02,plate002_F02,plate002,F02,i7-1789,TGGTCTGTTA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G02,plate002_G02,plate002,G02,i7-1758,TGCCAACATG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H02,plate002_H02,plate002,H02,i7-1643,TCATTACACG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A03,plate002_A03,plate002,A03,i7-0475,ATAACCTGAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B03,plate002_B03,plate002,B03,i7-1768,TGCGGTTCCA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C03,plate002_C03,plate002,C03,i7-0127,AATACTTGCC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D03,plate002_D03,plate002,D03,i7-0139,AATCGCGGAA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E03,plate002_E03,plate002,E03,i7-0661,CAGAACGCGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F03,plate002_F03,plate002,F03,i7-1722,TGAACCAAGG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G03,plate002_G03,plate002,G03,i7-1819,TGTTAGTCAG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H03,plate002_H03,plate002,H03,i7-0683,CAGTAGGTAA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A04,plate002_A04,plate002,A04,i7-0680,CAGGTACTTC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B04,plate002_B04,plate002,B04,i7-0678,CAGGAATATG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C04,plate002_C04,plate002,C04,i7-0043,AACGCACAAT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D04,plate002_D04,plate002,D04,i7-0576,ATTGAGAAGG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E04,plate002_E04,plate002,E04,i7-0578,ATTGCACCTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F04,plate002_F04,plate002,F04,i7-0116,AAGTGGATAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G04,plate002_G04,plate002,G04,i7-1192,GATAAGATGC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H04,plate002_H04,plate002,H04,i7-1555,TACCTCGACA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A05,plate002_A05,plate002,A05,i7-0028,AACCGAGCCA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B05,plate002_B05,plate002,B05,i7-1782,TGGATTCAAG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C05,plate002_C05,plate002,C05,i7-0586,CAACAGATAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D05,plate002_D05,plate002,D05,i7-0102,AAGGTAACTC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E05,plate002_E05,plate002,E05,i7-1728,TGACAATACG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F05,plate002_F05,plate002,F05,i7-1221,GATTGTGCAT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G05,plate002_G05,plate002,G05,i7-1231,GCACACCATT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H05,plate002_H05,plate002,H05,i7-1695,TCGTTATTCC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A06,plate002_A06,plate002,A06,i7-0597,CAACGTCATT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B06,plate002_B06,plate002,B06,i7-0068,AAGAACGATG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C06,plate002_C06,plate002,C06,i7-1059,CTGCAATTAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D06,plate002_D06,plate002,D06,i7-0907,CGCCGATGAT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E06,plate002_E06,plate002,E06,i7-0928,CGGACAAGAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F06,plate002_F06,plate002,F06,i7-0952,CGTACTCCTC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G06,plate002_G06,plate002,G06,i7-1757,TGCATGAGTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H06,plate002_H06,plate002,H06,i7-0141,AATCTGGAGC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A07,plate002_A07,plate002,A07,i7-0581,ATTGGTCAGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B07,plate002_B07,plate002,B07,i7-0075,AAGACCTGTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C07,plate002_C07,plate002,C07,i7-0917,CGCTAATGAA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D07,plate002_D07,plate002,D07,i7-0668,CAGAGTGCAT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E07,plate002_E07,plate002,E07,i7-1364,GGTAAGCTGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F07,plate002_F07,plate002,F07,i7-1432,GTATTCAGTG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G07,plate002_G07,plate002,G07,i7-1286,GCGCTACGTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H07,plate002_H07,plate002,H07,i7-1900,TTGTCATAGC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A08,plate002_A08,plate002,A08,i7-0558,ATGTTGCGGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B08,plate002_B08,plate002,B08,i7-0147,AATGCTAACC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C08,plate002_C08,plate002,C08,i7-1404,GTAACACGTA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D08,plate002_D08,plate002,D08,i7-1199,GATATACGGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E08,plate002_E08,plate002,E08,i7-0513,ATCCGAGAGG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F08,plate002_F08,plate002,F08,i7-0623,CAATTCACAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G08,plate002_G08,plate002,G08,i7-1899,TTGTCAGTTC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H08,plate002_H08,plate002,H08,i7-0485,ATACGCCATT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A09,plate002_A09,plate002,A09,i7-0636,CACCAATAAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B09,plate002_B09,plate002,B09,i7-1296,GCGTCCACAA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C09,plate002_C09,plate002,C09,i7-0555,ATGTGCGCTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D09,plate002_D09,plate002,D09,i7-1609,TATCTAGTGC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E09,plate002_E09,plate002,E09,i7-1268,GCCTACAATG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F09,plate002_F09,plate002,F09,i7-0096,AAGCTCAGTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G09,plate002_G09,plate002,G09,i7-1562,TACGCTTAGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H09,plate002_H09,plate002,H09,i7-0510,ATCCACTAGG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A10,plate002_A10,plate002,A10,i7-1802,TGTCCGTCTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B10,plate002_B10,plate002,B10,i7-0716,CATGAGTAAC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C10,plate002_C10,plate002,C10,i7-1063,CTGCGCGAAT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D10,plate002_D10,plate002,D10,i7-1074,CTGTAGTATG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E10,plate002_E10,plate002,E10,i7-1173,GAGCCGTACA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F10,plate002_F10,plate002,F10,i7-0727,CATTCTTAGG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G10,plate002_G10,plate002,G10,i7-1730,TGACCGACAA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H10,plate002_H10,plate002,H10,i7-1816,TGTGTAACCG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A11,plate002_A11,plate002,A11,i7-0862,CGAAGGACTG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B11,plate002_B11,plate002,B11,i7-1696,TCTACCGTCA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C11,plate002_C11,plate002,C11,i7-1365,GGTAATATCG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D11,plate002_D11,plate002,D11,i7-1813,TGTGCGAGTT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E11,plate002_E11,plate002,E11,i7-1848,TTATCGCTGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F11,plate002_F11,plate002,F11,i7-1400,GGTTGAGTTC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G11,plate002_G11,plate002,G11,i7-0580,ATTGGACGCC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H11,plate002_H11,plate002,H11,i7-1616,TATGTGTGTG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_A12,plate002_A12,plate002,A12,i7-1598,TAGTTATCGC,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_B12,plate002_B12,plate002,B12,i7-1755,TGCAGGTGAT,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_C12,plate002_C12,plate002,C12,i7-1809,TGTGAAGCTA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_D12,plate002_D12,plate002,D12,i7-0599,CAACTCCTGA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_E12,plate002_E12,plate002,E12,i7-1795,TGTAGCAACG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_F12,plate002_F12,plate002,F12,i7-0937,CGGTAACGCA,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_G12,plate002_G12,plate002,G12,i7-1423,GTAGCAGCAG,i5-0808,ATCTCCACGG,,
OCTOPUS_plate002_H12,plate002_H12,plate002,H12,i7-0989,CTACAGCCGA,i5-0808,ATCTCCACGG,,
```
</details>
<br>

# Day 2 - Run Wrapup

## 4. Glycerol Stock Culture Plates

1. Examine each Culture Plate for any evaporated wells or wells with little or no growth.
2. To your Culture Plate, add 100 µl 50% Glycerol to each well, mix, foil, and store at -80˚C.
3. If reusing lids from the plates, wash with 70% ethanol.

---
Notes:
- We use the Rainin Liquidator for adding glycerol and use the sterile lid from the tip box as a reservoir.
---

## 5. Run Analysis Pipeline

We provide an analysis pipeline for this workflow! Feel free to adapt it to your needs.

See the [Running the Analysis Pipeline](/docs/Running-the-Analysis-Pipeline.md) page for details.

