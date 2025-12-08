
# Sample trimming and variant calling analysis for Streptococcus pneumoniae in vivo evolution populations

# Trim reads
for i in $(cat ~/rawreads/sample_list.txt ) ; do mkdir ~/trimmed/"$i" ; done
for i in $(cat ~/rawreads/sample_list.txt ) ; do /opt/trimmomatic/Trimmomatic-0.36/bin/trimmomatic PE -threads 8 -phred33 -trimlog ~/trimmed/"$i"/trim.log ~/rawreads/"$i"/*R1_001.fastq.gz ~/rawreads/"$i"/*R2_001.fastq.gz ~/trimmed/"$i"/"$i".forward.trimmed.paired.fastq.gz ~/trimmed/"$i"/"$i".forward.trimmed.unpaired.fastq.gz ~/trimmed/"$i"/"$i".reverse.trimmed.paired.fastq.gz ~/trimmed/"$i"/"$i".reverse.trimmed.unpaired.fastq.gz ILLUMINACLIP:~/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee ~/trimmed/"$i"/"$i".trim.stdout ; done

# Run FastQC as a quality check
for i in $(cat ~/rawreads/sample_list.txt ) ; do mkdir ~/fastqc/"$i" ; /opt/fastqc/fastqc-0.11.5/fastqc -o ~/fastqc/"$i" -f fastq ~/trimmed/"$i"/"$i".forward.trimmed.paired.fastq.gz ~/trimmed/"$i"/"$i".forward.trimmed.unpaired.fastq.gz ~/trimmed/"$i"/"$i".reverse.trimmed.paired.fastq.gz ~/trimmed/"$i"/"$i".reverse.trimmed.unpaired.fastq.gz ; done

# breseq
for i in $(cat ~/rawreads/sample_list.txt ) ; do breseq -r ~/breseq/TIGR4vApr2021.gbk ~/trimmed/"$i"/"$i".forward.trimmed.paired.fastq.gz ~/trimmed/"$i"/"$i".reverse.trimmed.paired.fastq.gz ~/trimmed/"$i"/"$i".forward.trimmed.unpaired.fastq.gz ~/trimmed/"$i"/"$i".reverse.trimmed.unpaired.fastq.gz -o ~/breseq/"$i" -p -j 4 --polymorphism-minimum-variant-coverage-each-strand 3 --consensus-minimum-variant-coverage-each-strand 3 ; done


# Breseq output files from individual samples were combined into a single .xlsx file 
# Create a conda environment using Python 2.7
conda create -n breseq_parser_env python=2.7 -y
conda activate breseq_parser_env

# Install the required packages
conda install -c conda-forge \
    beautifulsoup4 \
    lxml \
    openpyxl \
    pandas \
    -y

python2 ~/scripts/BreseqCatEdited.py -p -d ~/breseq

cd ~/breseq

# Output files were covered to csvs
xlsx2csv Breseq_Output.xlsx -s 3 > breseq_output_nje.csv
xlsx2csv Breseq_Output.xlsx -s 2 > breseq_output_mc.csv
xlsx2csv Breseq_Output.xlsx -s 1 > breseq_output_snps.csv

# And converted to no utf characters
iconv -c -f utf-8 -t ascii//TRANSLIT breseq_output_snps.csv > breseq_output_snps_noutf.csv
iconv -c -f utf-8 -t ascii//TRANSLIT breseq_output_mc.csv > breseq_output_mc_noutf.csv
iconv -c -f utf-8 -t ascii//TRANSLIT breseq_output_nje.csv > breseq_output_nje_noutf.csv


#----- Get average coverage data for each sample ------------

# Copy output.gd file from every breseq output to a new directory and rename the output.gd file to "samplename.gd"
mkdir ~/breseq/outputfiles
for i in $(cat ~/rawreads/sample_list.txt) ; do cp ~/breseq/"$i"/output/output.gd ~/breseq/outputfiles/"$i".gd ; done

# navigate to outfiles directory
cd breseq/outputfiles

# get number of input bases (line 12 of output.gd file), becomes column 2
awk -v OFS=',' 'FNR==12{print FILENAME,$2}' *.gd > coverage_results.csv

# get number of mapped bases (line 14 of output.gd file), becomes column 3
awk -v OFS=',' 'FNR==14{print $2}' *.gd > coverage_results2.csv
paste --delimiters=',' coverage_results.csv coverage_results2.csv > coverage_results3.csv

# add column that contains genome length, becomes column 4
awk -F ',' -v OFS=',' '{ $(NF+1) = 2160842; print }' coverage_results3.csv > coverage_results4.csv

# add column that divides mapped bases by genome length
awk  -F "," -v OFS=',' '{$5 = $3 / 2160842}1' coverage_results4.csv > coverage_results5.csv

# add column that divides mapped bases by input bases
awk  -F "," -v OFS=',' '{$6 = $3 / $2}1' coverage_results5.csv > coverage_results6.csv

# add column names to file
echo "sample_name,input_bases,mapped_bases,genome_length,average_cov,percent_mapped_bases" > tmp
cat coverage_results6.csv >> tmp
mv tmp coverage_results.csv 

# remove intermediate files
rm coverage_results2.csv
rm coverage_results3.csv
rm coverage_results4.csv
rm coverage_results5.csv
rm coverage_results6.csv
