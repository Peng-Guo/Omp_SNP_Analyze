---
title: "Omp Data Analysis Report"
author: "guopeng"
data: 20240313
output:
  html_document:
    df_print: paged
  rmdformats::robobook:
    highlight: tango
---

# 数据获取

## 数据筛选

### Pipeline

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
grep "Escherichia coli" assembly_summary.txt > ecoli.all.txt
wc -l ecoli.all.txt
> 260227 ecoli.all.txt

awk -F'\t' '$21 == "na" {print $1 "\t" $20}' ecoli.all.txt > ecoli.refseq.txt
wc -l ecoli.refseq.txt
> 37360 ecoli.refseq.txt

head -2 ecoli.refseq.txt
> GCA_002442165.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/442/165/GCA_002442165.1_ASM244216v1
> GCA_002442365.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/442/365/GCA_002442365.1_ASM244236v1
```

## 数据下载

### Pipeline

```bash
awk -F'\t' '{print $2}' ecoli.refseq.txt > ecoli.ftp.list
mkdir ecoli_genome_gz
cd ecoli_genome_gz
bash ../download.ftp.sh ../ecoli.ftp.list
ls > ../ecoli.genome.gz.list
```

### Script

- download.ftp.sh

```bash
for link in `cat $1`
do
	name=$(echo $link | awk -F'/' '{print $NF "_genomic.fna.gz"}')
	wget "${link}/${name}"
done
```

# 基因预测

## 数据解压

### Pipeline

```bash
mkdir ../ecoli_genome_fasta
bash ../ungz.sh ../ecoli.genome.gz.list ../ecoli_genome_fasta
cd ../ecoli_genome_fasta
ls > ../ecoli.genome.fasta.list
```

### Script

- ungz.sh

```bash
for file in `cat $1`
do
	path="$2/$(basename ${file} .fna.gz).fasta"
	gunzip -c ${file} > ${path}
done
```

## 基因预测

### Pipeline

```bash
mkdir ../ecoli_gene_predict
bash ../gene.predict.sh ../ecoli.genome.fasta.list ../ecoli_gene_predict
```

### Script

- prodigal.sh

```bash
mkdir $2/gff
mkdir $2/gene
mkdir $2/protein

for file in `cat $1`
do
	filename=${file%.*}
	prodigal -i $file -o "$2/gff/${filename}.gff" -d "$2/gene/${filename}.gene.fasta" -a "$2/protein/${filename}.protein.fasta" -p single
done
```

# blastp

## 基因合并

### Pipeline

```bash
cd ../ecoli_gene_predict/protein
ls > ../ecoli.gene.protein.list
cd ..
bash ../combined.sh ecoli.gene.protein.list protein > ecoli.gene.protein.summary.fasta

cd gene
ls > ../ecoli.gene.list
cd ..
bash ../combined.sh ecoli.gene.protein.list protein
```

### Script

- combined.sh

```bash
for file in `cat $1`
do
	cat "$2/${file}"
done
```

## 下载参考序列

### Pipeline

```bash
mkdir ../omp_prot_db
cd ../omp_prot_db
wget -O ompa.prot.ref https://rest.uniprot.org/uniprotkb/P0A910.fasta
wget -O ompc.prot.ref https://rest.uniprot.org/uniprotkb/P06996.fasta
wget -O ompf.prot.ref https://rest.uniprot.org/uniprotkb/P02931.fasta
cat omp*.prot.ref >> omp.prot.ref.fasta
```

## 建库与运行

### Pipeline

```bash
diamond makedb --in omp.prot.ref.fasta --db omp.prot
diamond blastp -d omp.prot -q ../ecoli_gene_predict/ecoli.gene.protein.summary.fasta -o blastp.out
```

# 序列筛选

## 基因分类

### Pipeline

```bash
mkdir ../omp_prot_fasta
cd ../omp_prot_fasta
bash ../gene.filter.sh omp.gene.list ../ecoli_gene_predict/ecoli.gene.protein.summary.fasta > omp.prot.fasta
```

[Get `omp.gene.list`](./omp.filter.python.html)

### Script

- gene.filter.sh

```bash
for name in `cat $1`
do
        awk "/${name} #/,/\\*/" "$2"
done
```

## 基因去冗

### Pipeline

```bash
cd-hit -i omp.prot.fasta -o omp.uniq -c 1.0 -n 5 -d 0 -M 150000 -T 52 -aS 1.0 -aL 1.0 -g 1
python omp.prot.freq.py omp.uniq.clstr omp.prot.freq
grep -c '^1\s' > omp.prot.filted
grep -Po ', \K[^at]*' omp.prot.filted | sed 's/....$//' > omp.uniq.filted.list
bash ../gene.filter.sh omp.uniq.filted.list ../ecoli_gene_predict/ecoli.gene.summary.fasta > omp.gene.fasta
```

### Script

```python
import pandas as pd
import numpy as np
import sys, os

if __name__ == '__main__':
    infile = sys.argv[1] #
    outfile = sys.argv[2]


    repre_seqs = {}
    seqids = {}
    seqsamples = {}
    with open(infile, 'r') as fin:
        for l in fin.readlines():
            seqflag = 0
            if l.startswith('>'):
                cluster = l.strip('>').strip('\n').split(" ")[1]
            else:
                seqid = l.split('>')[1].split('|')[-1].split('.')[0]
                sampleid = seqid.split('_')[0]
                if len(seqids) == 0:
                    seqids[cluster] = seqid
                else:
                    if cluster in seqids.keys():
                        seqids[cluster] = seqids[cluster] + ',' + seqid
                    else:
                        seqids[cluster] = seqid
                if "*" in l:
                    repre_seqs[cluster] = seqid

    outdf = pd.DataFrame({'Cluster':list(seqids.keys()), 'represent_seqid':list(repre_seqs.values()), 'seqid':list(seqids.values())})

    #outdf = pd.merge(outdf, summarydf, on = 'represent_seqid', how = 'left')

    outdf.to_csv(outfile, sep = '\t', index=False)
```

# SNP分析

## 生成VCF

### Pipeline

```bash
mkdir ../omp_snp
cd ../omp_snp
bash snp.sh ompa.gene.ref ompa.gene.fasta ./ompa/ompa
```

### Script

- snp.sh

```bash
minimap2 -ax asm20 $1 $2 > "$3.sam"
samtools view -S -b "$3.sam" > "$3.bam"
samtools sort "$3.bam" -o "$3.sorted.bam"
bcftools mpileup -Ou -f $1 "$3.sorted.bam" | bcftools call -mv -Ov -o "$3.vcf"
```


## 建库与注释

### Pipeline

```bash
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip


cd snpEff  #进入 snpEff 目录下
mkdir data  #新建 data 目录
cd data  #进入 data 目录下，必须在该目录
mkdir genomes  #新建 genomes 目录，用于建立 ecoli 目录中的 bin
mkdir ecoli  #新建 ecoli 目录,对应物种名或后续软件调用的参数名

#将参考基因序列放至 genomes 目录下并改名，改名对应 data 中的 ecoli 文件夹名
gzip -d {reference}_genomic.fna.gz
mv dir_path/{reference}_genomic.fna snpEff/data/ecoli.fa

#将参考注释放至 ecoli 文件夹中并一定改名为 gene.gff
gzip -d {reference}_genomic.gff.gz
mv dir_path/{reference}_genomic.gff snpEff/data/genes.gff

#配置 snpEff/snpEff.config 文件
echo "ecoli.genome:ecoli" >> snpEff.conf
#可配置 config 文件自定义 data 文件夹路径
#data_dir = PATH/TO/data/

#软件配置自动建立数据库，建立后文件夹 ecoli 中生成 snpEffectPredictor.bin 文件
java -jar snpEff/snpEff.jar build -gff3 ecoli
#-gff3 对应为 ecoli 文件夹中注释文件的对应格式版本
#如果注释文件的格式版本在后缀名中无法体现，可通过 head 命令查看文件头信息获得
#head ecoli.gff
#ecoli 对应为 data 目录中的 ecoli 文件夹以及 genomes 目录中的 ecoli.fa 文件

#配置完成后的目录结构，其他文件夹/结构不显示
snpEff
├── SnpSift.jar
├── data
│   ├── ecoli
│   │   ├── genes.gff
│   │   └── snpEffectPredictor.bin
│   ├── genomes
│   │   └── ecoli.fa
├── examples/
├── galaxy/
├── scripts/
├── snpEff.config
└── snpEff.jar
```

