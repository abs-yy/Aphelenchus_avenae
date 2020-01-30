# Introduction
For many reasons, I was looking for genomes of anhydrobiotic nematodes when I found sequencing data (DNA and RNA) of Aphelenchus avenae. There were genome/transcriptome assemblies and no papers assosiated with this study.

# Method
## Data Acquisition
- Sequence data list
  - DNA-Seq data
    - 
  - RNA-Seq data
    - SRR1175676, SRR1175692, SRR1175695, SRR1175696, SRR1175697, SRR1175706, SRR1175707, SRR1175708, SRR1175729, SRR1175731, SRR1175736, SRR1175737, SRR1175739, SRR1175740, SRR1174913
- Data Acquisition
  - sratools prefetch
  - sratools fastq-dump
  - Data filtering
    - There were reads that pairs could not be found in several DNA-Seq data, so I used fastq-pair to filter the data
      -  https://github.com/linsalrob/fastq-pair/blob/master/CITATION.md
      - /path/to/fastq-pair/build/fastq_pair SRX476031.sra_2.fastq SRX476031.sra_4.fastq
  
## Genome assembly
- There are many genome assemblers available, I used two famous softwares to test if there are any differences.
- SPADES v3.14.0
  ```
  ~/bin/SPAdes-3.14.0-Linux/bin/spades.py \
    --pe-1 1 SRR1180010.1.sra_1.fastq \
    --pe-2 1 SRR1180010.1.sra_2.fastq \
    --pe-1 2 SRR1176881.sra_1.fastq.paired.fq \
    --pe-2 2 SRR1176881.sra_2.fastq.paired.fq \
    --pe-1 3 SRR1179837.sra_1.fastq \
    --pe-2 3 SRR1179837.sra_2.fastq \
    --mp-1 2 SRX476031.sra_2.fastq.paired.fq \
    --mp-2 2 SRX476031.sra_4.fastq.paired.fq \
    --mp-1 3 SRR1176816.sra_2.fastq.paired.fq \
    --mp-2 3 SRR1176816.sra_4.fastq.paired.fq \
    -o Aphelenchus_avenae_spades -t 64
  ```
- MaSuRCA v3.3.5
  - Commands
    ```
    % /path/to/MaSuRCA-3.3.5/bin/masurca -g 
    # edit configureation file as shown below
    % /path/to/MaSuRCA-3.3.5/bin/masurca config.txt
    % ./assemble.sh
    ```
  - Config.file
  ```
    # example configuration file

    # DATA is specified as type {PE,JUMP,OTHER,PACBIO} and 5 fields:
    # 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads
    # 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be
    # innies, i.e. --->.<---, and JUMP are assumed to be outties
    # <---.--->. If there are any jump libraries that are innies, such as
    # longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads
    # are optional for PE libraries and mandatory for JUMP libraries. Any
    # OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
    # converted into Celera Assembler compatible .frg files (see
    # http://wgs-assembler.sourceforge.com)
    DATA
    #Illumina paired end reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
    #if single-end, do not specify <reverse_reads>
    #MUST HAVE Illumina paired end reads to use MaSuRCA
    PE= p1 500 50 /path/to/data/SRR1180010.1.sra_1.fastq /path/to/data/SRR1180010.1.sra_2.fastq
    PE= p2 500 50 /path/to/data/SRR1176881.sra_1.fastq.paired.fq /path/to/data/SRR1176881.sra_2.fastq.paired.fq
    PE= p3 500 50 /path/to/data/SRR1179837.sra_1.fastq /path/to/data/SRR1179837.sra_2.fastq
    #Illumina mate pair reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
    #JUMP= sh 3600 200  /FULL_PATH/short_1.fastq.paired.fq  /FULL_PATH/short_2.fastq.paired.fq
    JUMP= j1 20000 1000 /path/to/data/SRX476031.sra_2.fastq.paired.fq /path/to/data/SRX476031.sra_4.fastq.paired.fq
    JUMP= j2 8000 500 /path/to/data/SRR1176816.sra_2.fastq.paired.fq /path/to/data/SRR1176816.sra_4.fastq.paired.fq
    #pacbio OR nanopore reads must be in a single fasta or fastq file with absolute path, can be gzipped
    #if you have both types of reads supply them both as NANOPORE type
    #PACBIO=/FULL_PATH/pacbio.fa
    #NANOPORE=/FULL_PATH/nanopore.fa
    #Other reads (Sanger, 454, etc) one frg file, concatenate your frg files into one if you have many
    #OTHER=/FULL_PATH/file.frg
    #synteny-assisted assembly, concatenate all reference genomes into one reference.fa; works for Illumina-only data
    #REFERENCE=/FULL_PATH/nanopore.fa
    END

    PARAMETERS
    #PLEASE READ all comments to essential parameters below, and set the parameters according to your project
    #set this to 1 if your Illumina jumping library reads are shorter than 100bp
    EXTEND_JUMP_READS=0
    #this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
    GRAPH_KMER_SIZE = auto
    #set this to 1 for all Illumina-only assemblies
    #set this to 0 if you have more than 15x coverage by long reads (Pacbio or Nanopore) or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc)
    USE_LINKING_MATES = 0
    #specifies whether to run the assembly on the grid
    #USE_GRID=1
    #specifies grid engine to use SGE or SLURM
    #GRID_ENGINE=SGE
    #specifies queue (for SGE) or partition (for SLURM) to use when running on the grid MANDATORY
    #GRID_QUEUE=all.q
    #batch size in the amount of long read sequence for each batch on the grid
    #GRID_BATCH_SIZE=500000000
    #use at most this much coverage by the longest Pacbio or Nanopore reads, discard the rest of the reads
    #can increase this to 30 or 35 if your reads are short (N50<7000bp)
    LHE_COVERAGE=25
    #set to 0 (default) to do two passes of mega-reads for slower, but higher quality assembly, otherwise set to 1
    MEGA_READS_ONE_PASS=0
    #this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms
    LIMIT_JUMP_COVERAGE = 300
    #these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically.
    #CABOG ASSEMBLY ONLY: set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
    CA_PARAMETERS =  cgwErrorRate=0.15
    #CABOG ASSEMBLY ONLY: whether to attempt to close gaps in scaffolds with Illumina  or long read data
    CLOSE_GAPS=1
    #auto-detected number of cpus to use, set this to the number of CPUs/threads per node you will be using
    NUM_THREADS = 64
    #this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20
    #JF_SIZE = 200000000
    JF_SIZE=6644521840
    #ILLUMINA ONLY. Set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>=8Gbp) genomes from Illumina-only data
    SOAP_ASSEMBLY=0
    #Hybrid Illumina paired end + Nanopore/PacBio assembly ONLY.  Set this to 1 to use Flye assembler for final assembly of corrected mega-reads.  A lot faster than CABOG, at the expense of some contiguity. Works well even when MEGA_READS_ONE_PASS is set to 1.  DO NOT use if you have less than 15x coverage by long reads.
    FLYE_ASSEMBLY=0
    END
    ```

## Transcriptome assembly
