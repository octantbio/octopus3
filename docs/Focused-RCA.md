**Focused RCA (optional)**
==========================

One modification made to the RCA process ([step 1.2](/docs/Bench-Protocol.md#12-rca-reaction)) was to contain, in addition to the random hexamers, "focused" primers specific to the plasmid samples that allow for enhanced amplification of plasmid DNA as opposed to gDNA or other contaminants, leading to more coverage of the relevant data. We observe that without the focused primers, genomic host DNA would represent ~50% of the reads. With the focused primers, genomic host DNA is reduced to under 10%. Implementing these primers in your RCA protocol is optional, though we find that it significantly reduces dropouts from samples with borderline coverage. We note that this is not always possible if your samples don't all contain a common sequence and that other options are available (e.g. sequencing deeper).

Below are the exact custom exo-resistant primers we ordered from [IDT](https://www.idtdna.com/pages/support/faqs/how-do-i-order-a-phosphorothioate-modified-oligo-) specific to plasmids that contain the ColE1 origin of replication. Loosely, we designed them to be ~15bp, with a Tm of ~40 °C, with ~40% GC content, ending in a G or C, and somewhat evenly spaced out across the common ColE1 sequence. The last three 3' bases are phosphorothiorated to make them resistant to phi29 3'->5' exonuclease activity.

### **Table of "focused" primers for ColE1 origin of replication**
Asterisk (*) between bases denote a phosphorothioate chemical modification.

|      Name     	| Sequence (5' -> 3')     	| Scale 	| Tm 	| Primer Len 	|
|:-------------:	|-------------------------	|:-----:	|:--:	|:----------:	|
| Fwd-ORI-RCA-1 	| AAGGATCTTCTT\*G\*A\*G   	| 100nm 	| 40 	|     15     	|
| Fwd-ORI-RCA-3 	| TCAAGAGCTACC\*A\*A\*C   	| 100nm 	| 44 	|     15     	|
| Fwd-ORI-RCA-4 	| CAAATACTGTTC\*T\*T\*C   	| 100nm 	| 37 	|     15     	|
| Fwd-ORI-RCA-2 	| GTCGGAACAGGA\*G\*A\*G   	| 100nm 	| 48 	|     15     	|
| Rev-ORI-RCA-1 	| CTGGTAACAGGA\*T\*T\*A   	| 100nm 	| 40 	|     15     	|
| Rev-ORI-RCA-2 	| AAGACACGACTT\*A\*T\*C   	| 100nm 	| 41 	|     15     	|
| Rev-ORI-RCA-3 	| TATCCGGTAACT\*A\*T\*C   	| 100nm 	| 39 	|     15     	|
| Rev-ORI-RCA-4 	| GACTATAAAGATAC\*C\*A\*G 	| 100nm 	| 39 	|     17     	|

<details>
<summary>copy-paste friendly version of sequences</summary>

```
AAGGATCTTCTT*G*A*G
TCAAGAGCTACC*A*A*C
CAAATACTGTTC*T*T*C
GTCGGAACAGGA*G*A*G
CTGGTAACAGGA*T*T*A
AAGACACGACTT*A*T*C
TATCCGGTAACT*A*T*C
GACTATAAAGATAC*C*A*G
```
</details>
<br>

### **Template sequence for the ColE1 origin of replication**
```
>Ori_ColE1 [ident:99.84, cov:100.00]
ATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCACCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAA
```

![ColE1 origin of replication with "focused" primers](/img/cole1-ori-focused-primers.png)

