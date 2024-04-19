# File Tree

- assembly_summary.txt
- ecoli.ftp.link
- download.ftp.sh
- ecoli.gene.list
- prodigal.sh
- **ecoli_gene_gz**
- **prodigal_gff**
- **prodigal_gene**
- **prodigal_protein**

# Pipeline

1. 通过FTP下载e.coli基因组数据

下载`assembly_summary.txt`，其中包含了genbank中bacteria的所有assembly数据；

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
```

从`assembly_summary.txt`筛选*E.coli*的FTP链接；

```bash
awk -F'\t' '$8 == "Escherichia coli" {print $20}' assembly_summary.txt > ecoli.ftp.link
```

根据`ecoli.ftp.link`文件下载*E.coli*的基因组数据；

```bash
mkdir ecoli_gene_gz
cd ecoli_gene_gz
bash ../download.ftp.sh ../ecoli.ftp.link
```

> download.ftp.sh内容如下：
>
> ```bash
> for link in `cat $1`
> do
> 	name=$(echo $link | awk -F'/' '{print $NF "_genomic.fna.gz"}')
> 	wget "${link}/${name}"
> done
> ```
>

2. 使用prodigal预测基因组上的基因，得到基因的核苷酸序列和氨基酸序列

解压得到所有的基因组fasta文件；

```bash
bash ungz.sh ecoli_gene_gz ecoli_gene_fasta
```

获取下载得到的基因名称列表；

```bash
cd ecoli_gene_fasta
ls > ../ecoli.gene.list
```

使用prodigal进行基因预测；

```bash
mkdir ecoli_gene_fasta
bash prodigal.sh ecoli.gene.list ecoli_gene_fasta
```

> ungz.bash内容如下：
>
> ```bash
> for file in "$1"/*;
> do
> 	path="$2/$(basename "$file" .fna.gz).fasta"
>   gunzip -c "$file" > "$path"
> done
> ```
>
> prodigal.sh内容如下：
>
> ```bash 
> mkdir prodigal_gff
> mkdir prodigal_gene
> mkdir prodigal_protein
> 
> for file in `cat $1`
> do
>   filename=${file%.*}
>   prodigal -i "$2/$file" -o "./prodigal_gff/${filename}.gff" -d "./prodigal_gene/${filename}.gene.fasta" -a "./prodigal_protein/${filename}.protein.fasta" -p single
> done
> ```
>

3. 合并所有基因的氨基酸序列

获取下载得到的氨基酸序列名称列表；

```bash
cd prodigal_protein
ls > ../ecoli.protein.list
```

合并所有的氨基酸序列为一个fasta文件；

```bash
bash combined.sh ecoli.protein.list prodigal_protein
```

> combined.sh内容如下：
>
> ```bash 
> for file in `cat $1`
> do
> 	cat "$2/${file}" >> all.protein.fasta
> done
> ```
>

4. 下载omp基因的参考核苷酸和氨基酸序列

```txt
ompa: https://www.uniprot.org/uniprotkb/P0A910/entry
ompc: https://www.uniprot.org/uniprotkb/P06996/entry
ompf: https://www.uniprot.org/uniprotkb/P02931/entry
```



5 使用omp基因的氨基酸序列进行Diamond blast建库（mmseq）



6 对合并得到的所有氨基酸序列进行Diamond blastp



7 根据Diamond结果对所有基因进行初步的分类和筛选



8 对筛选结果使用cd-hit去除重复序列，得到unique氨基酸序列



9 根据unique氨基酸序列筛选对应的unique核苷酸序列



10 使用minimap2以参考核苷酸序列为基础，对unique核苷酸序列进行比对，得到sam文件



11 使用samtools对sam文件进行转化，排序，得到bam文件



12 使用bcftools处理bam文件，得到vcf文件



13 使用snpeff对vcf文件进行注释，展示snp位点，得到anno.vcf文件
