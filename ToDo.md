

# TODO: Issues and Fixes for humann4_tools

## 🐞 Bug Fixes
- [ ] move kneaddata and humann4 output files to a different directory then delete the temp directory
- [ ] in cli.py , if args.use_metadata is True but no --seq-dir is provided, the current code does not execute metadata-based file collection.
- [ ] add the flags and examples for skiip kneaddata	
- [ ] the logic for selecting kneaddata files isn't working. pulls in all samples

`

## 🔧 Features to Improve
- [ ] Check if intermediate files already exist before performing steps in the analysis
        e.g. if Kneaddata output files already exist, skip kneaddata
- [ ] Add genedirectory, pathwaydirectory, etc flags to options section in README.md
- [ ] incorporate into the README.md new join_unstratify_humann_output command


## 📅 Future Plans
- [ ] Test SHAP force plots with new dataset


```bash


	
# not parallel	
	humann4-tools --run-preprocessing \
	--log-file test.log \
	--log-level DEBUG \
	--output-dir test.out2 \
	--output-prefix test- \
	--kneaddata-dbs  /home/david/Databases/KneadDataDB /home/david/Databases/BT2ContaminantDB \
	--threads 24 \
	--sample-key TestKey.csv \
	--use-metadata \
	--seq-dir /home/david/Data/TrimmedMSSFiles/ \
  	--kneaddata-output-dir ~/Data/TrimmedMSSFiles \
  	--humann4-output-dir ~/Documents/Alignments/Humann4Alignments \
	--pathway-dir test.out2 \
	--gene-dir test.out2 \
	--group-col "Group" \
	--paired \
	--decontaminate-pairs strict \
	--annotations-dir annotations 


# parallel
humann4-tools 	--run-preprocessing \
	--log-file test.log \
	--log-level DEBUG \
	--output-dir test.out2 \
	--output-prefix test- \
	--kneaddata-dbs  /home/david/Databases/KneadDataDB /home/david/Databases/BT2ContaminantDB \
	--threads 24 \
	--sample-key TestKey.csv \
	--use-metadata \
	--seq-dir /home/david/Data/TrimmedMSSFiles/ \
  	--kneaddata-output-dir ~/Data/TrimmedMSSFiles \
  	--humann4-output-dir ~/Documents/Alignments/Humann4Alignments \
	--pathway-dir test.out2 \
	--gene-dir test.out2 \
	--group-col "Group" \
	--paired \
	--use-parallel \
	--threads-per-sample 4 \
	--decontaminate-pairs strict \
	--annotations-dir annotations 

# join, normalize, split stratify only
humann4-tools --join-only \
  --sample-key TestKey.csv \
  --pathway-dir /path/to/pathways \
  --gene-dir /path/to/genes \
  --output-dir joined_output \
  --units cpm \
  --no-interactive

  # skip kneaddata

  # Basic usage with skip-kneaddata and automatic file detection
humann4-tools --run-preprocessing --skip-kneaddata --input-fastq sample1.fastq sample2.fastq

# Using specific KneadData output files
humann4-tools --run-preprocessing --skip-kneaddata --kneaddata-output-files /path/to/sample1_paired_1.fastq /path/to/sample1_paired_2.fastq 

# Using a pattern to find KneadData output files
humann4-tools --run-preprocessing --skip-kneaddata --kneaddata-output-pattern "/path/to/kneaddata_output/*/*paired*.fastq"

# Using a pattern with sample names
humann4-tools --run-preprocessing --skip-kneaddata --kneaddata-output-pattern "/path/to/kneaddata_output/{sample}*paired*.fastq"
