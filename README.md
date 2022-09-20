# 1000MonkeyGenome
This repository contains scripts related to manuscript "Adaptive Structural Variants Contribute to Human Brain Development Revealed by 1,026 Rhesus Macaque Genomes"

# Table of contents
- Content
- Contact

## Content

- **process_maf.py:** An inhouse script was developed to process the Multiple Alignment Format (MAF) file and get the details of alignments of inversions and flanking regions.
      Method detail: Ten NHPs were retained for lineage tracing for inversions between humans and rhesus macaques, including chimpanzees (panTro6), gorillas (gorGor6), orangutans (ponAbe3), gibbons (nomLeu3), crab-eating macaques (macFas5), baboons (papAnu4), green monkeys (chlSab2), marmosets (calJac4), squirrel monkeys (saiBol1) and rhesus macaques (rheMac10Plus). We first downloaded the Multiple Alignment Format (MAF) file derived from syntenic net alignment for each NHP from the UCSC Genome Browser except for crab-eating macaques and rhesus macaques. For crab-eating macaques and rhesus macaques, the MAF file was not available in the UCSC Genome Browser, so we manually generated it using netToAxt (https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netToAxt). 
      
      # Criteria used to determine inversion alignments
      
      (1) use flanking regions upstream and downstream of inversion, which are twice the length of inversion as background
      
      (2) use 50% as a cutoff of fraction when determining chromosome
      
      (3) use 80% as a cutoff of fraction when determining alignment direction

      
      # Usage (Python) : python get_invMaf_aln.py \
      --inv hrInv.hg38.bed \
      --maf hg38.$species.synNet.maf.gz --flankFrac 2 \
      --chrSize $chrSize \
      --intersect inv.body.maf.hg38_${species}.list  \
      --intersect1 inv.up.maf.hg38_${species}.list \
      --intersect2 inv.down.maf.hg38_${species}.list \
      --processMaf processedMaf.hg38_${species}.list
  
      Options:
    
        Input:
        --inv: Input file (Inversion list. BED format in hg38 coordinates)
    
        --maf: MAF files between human and other primates
      
        --chrSize: Chromosome size file (Column 1： chr ID; column 2: chrmosome size) 
      
        Output:
        --intersect: intermediate file for alignments of inversion body region
      
        --intersect1: intermediate file for alignments of inversion upstream regions
      
        --intersect2: intermediate file for alignments of inversion downstream regions
      
        --processMaf: 
      
- **determine_inv.py:** Determine inversion status for the output of process_maf.py

      # Criteria used to determine inversion status
      
      (1) at least one flanking region fell on the same chromosome as the inversion body region 
      
      (2) these two regions were aligned in the opposite direction with respect to human
      
      (3) inversions were removed from downstream ancestor state reconstruction if they were not detected as inverted under criteria (1)(2) between human and rhesus macaque 


      # Usage (Python) : python determine_inv.py \
      --inv hrInv.hg38.bed \
      --maf inv.body.maf.hg38_${species}.list \
      --mafUp inv.up.maf.hg38_${species}.list \
      --mafDown inv.down.maf.hg38_${species}.list \
      --output invStatus.hg38_${species}.flankFrac2.noContext.txt
     
      Options:
    
        Input:
        --inv: Input file (Inversion list. BED format in hg38 coordinates)
    
        --maf: output file for alignments of inversion body region for script process_maf.py
      
        --mafUp:  output file for alignments of inversion upstream region for script process_maf.py
      
        --mafDown: output file for alignments of inversion downstream region for script process_maf.py
      
        Output:
        --output: inversion status for input file between humand and corresponding species


## Contact

- dingwq@pku.edu.cn
- 2101112223@pku.edu.cn
- zhangjie_imm@pku.edu.cn
