[TOC]
# 重测序与群体遗传
## 背景
- 基因组与比较基因组：
    - 基因组学（侠义：构建一个物种的基因组图谱）
    - 比较基因组学（物种之间的差异）
        - 物种之间层次的比较：那些物种近，哪些物种远，这些物种的分化时间
        - 基因组层次的比较：物种之间哪些基因发生了重组
        - 基因家族层次的比较：某一个物种的基因家族比另一个物种多了几条基因，某个基因在某个物种演化速率更快一些，正选择
- 重测序（已经有了参考基因组）：
    - 个体重测序：**本质是基因分型**（每个位点基因型是什么）（一般是应用）
    - 群体重测序（同一个物种的不同种类分别测序）：**本质是等位基因频率**
        - 研究进化关系
        - 基因定位：表型背后的基因位点
        - GWAS：全基因组关联分析
## 重测序分析一般步骤：
- 测序数据 => 数据质控 => 基因组比对 => (SNV,InDel,SNP,SV) => (群体进化分析，群体选择分析，GWAS分析)
- SNP单核苷酸多态性
- SNV单碱基变异
- CNV拷贝数变异
## 测序的应用
- 测序应用
    - 序列构建
        - DNA参考基因组构建
            - 物种参考基因组构建
            - 物种的Pan Genome
            - 宏基因组
        - RNA转录组构建
    - 变异检测
        - DNA
            - 按变异类型
                - SNP
                - INDEL
                - SV
            - 按测序范围
                - 全基因组WGS
                - 外显子组WES
                - 目标基因TGS
                - 酶切位点RAD/GBS
            - 按研究对象
                - 个体-不同个体之间的差异（遗传病）
                - 群体-群体等位基因频率的差异（群体遗传病）
                - 细胞-不同细胞之间的差异（肿瘤）
        - RNA
    - 丰度估计
        - DNA：按研究对象
            - 染色体：21三体综合征
            - 基因：拷贝数变异
        - RNA：表达丰度估计
# loading ===>
```
└── 01.variant_calling
    ├── 01.mapping
    ├── 02.SNP_indel
    ├── data
    └── ref
```
## step1_variant_calling
### mapping
- mapping
    - 构建index：构建输入参考基因组构建bwa，samtools，picard的index
    - bwa比对samtools排序：
    - picard去重：
```
$ cat ref/index.sh
#!/bin/bash

ln -s ../data/genome.fasta ./genome.fasta
samtools faidx genome.fasta
bwa index genome.fasta
picard CreateSequenceDictionary R=genome.fasta

$ cat 01.mapping/step1.bwa.sort.sh
#!/bin/bash

for i in `seq 2`;do
        bwa mem -t 2 -R "@RG\tID:S${i}\tSM:S${i}\tPL:illumia" ../ref/genome.fasta ../data/S${i}_1.fq.gz ../data/S${i}_2.fq.gz | samtools sort -@ 2 -m 1G -o S${i}.sort.bam -
done

# bwa两种比对模式，mem比较快
# -t 线程
# -R reads分组信息
# -@ 线程
# -m 内存
# -o 输出

$ cat 01.mapping/step2.picard_MarkDuplicates.sh
#!/bin/bash

# bwa比对后生成的sam文件用samtools sort排序并压缩成bam文件，用picard对排序好的bam文件去重
for i in `seq 2`;do
        picard -Xmx4g MarkDuplicates I=S${i}.sort.bam O=S${i}.sort.rmdup.bam CREATE_INDEX=true REMOVE_DUPLICATES=true M=S${i}.marked_dup_metrics.txt
done

# CREATE_INDEX=true创建bam文件的index，IGV,GATK等需要
# M=S${i}.marked_dup_metrics.txt 去重的日志信息
```
### SNP_InDel calling
- SNP calling
    - 哪些位点有单碱基多态性
    - 哪些位点相对于考基因组存在至少一个变异碱基（variant calling）
- genotype calling
    - 确定每个个体基因型的过程

- gatk：
    - 输入：之前比对好排完序去完重的bam文件，和参考基因组文件
    - 输出：.g.vcf文件（包含突变和未突变的位点）
```
$cat 02.SNP_indel/GATK_pipeline.sh
#!/bin/bash

# step1_HaplotypeCaller
for i in `seq 2`;do
        [ -d tmp ] || mkdir tmp
        gatk --java-options "-Xmx20g -Djava.io.tmpdir=./tmp" HaplotypeCaller -R ../ref/genome.fasta -I ../01.mapping/S${i}.sort.rmdup.bam -ERC GVCF -O S${i}.g.vcf &> S${i}.HC.log
        # -R 参考基因组 samtools创建的fail.index和picard的dict
        # -ERC GVCF 生成GVCF文件，如果不指定直接生成vcf文件
        # -O gvcf输出名称
done

# step2_combineGVCFs:并gvcf
# 将所有的gvcf生成一个表
ls ./S*.g.vcf > gvcf.list
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" CombineGVCFs -R ../ref/genome.fasta -V gvcf.list -O all.merge.g.vcf

# step3_genotypeGVCFs:每个个体Genotype的信息
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" GenotypeGVCFs -R ../ref/genome.fasta --variant all.merge.g.vcf -O all.merge_raw.vcf

# step4_filter_SNP:
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" SelectVariants -R ../ref/genome.fasta -V all.merge_raw.vcf --select-type SNP -O all.raw.snp.vcf
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" VariantFiltration -R ../ref/genome.fasta -V all.raw.snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter' -O all.filter.snp.vcf
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" SelectVariants -R ../ref/genome.fasta -V all.filter.snp.vcf --exclude-filtered -O all.filtered.snp.vcf

# step5_filter_indel:
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" SelectVariants -R ../ref/genome.fasta -V all.merge_raw.vcf --select-type INDEL -O all.raw.indel.vcf
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" VariantFiltration -R ../ref/genome.fasta -V all.raw.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'INDEL_filter' -O all.filter.indel.vcf
gatk --java-options "-Xmx10g -Djava.io.tmpdir=./tmp" SelectVariants -R ../ref/genome.fasta -V all.filter.snp.vcf --exclude-filtered -O all.filtered.indel.vcf
```
### SV calling
- 软件：lumpy
- 原理：read pair；split reads；read depth
- 输入：排完序去完重的bam文件
- 运行过程：
    - 提取非正常成对比对的reads：`samtools view -F 1294`
    - 提取split reads：lumpy中的`extractSplitReads_BwaMem`
    - 群体水平SV检测：`lumpyexpress`
    - 个体genotype calling：
        - 提取个体vcf`vcftools`
        - 进行检测：`svtype`
        - 合并个体vcf：`svtools`
```
$ cat SV_pipeline.sh
#!/bin/bash

lumpy="/pub/software/lumpy-sv"

for i in `seq 2`;do
        # step1_dicordants_bam
        samtools view -b -F 1294 -@ 12 ../01.mapping/S${i}.sort.rmdup.bam > S${i}.discordants.bam

        # step2_splitters_bam
        samtools view -h ../01.mapping/S${i}.sort.rmdup.bam | $lumpy/scripts/extractSplitReads_BwaMem -i stdin | samtools sort - > S${i}.splitters.bam
done


# step3_luppy
$lumpy/bin/lumpyexpress -B ../01.mapping/S1.sort.rmdup.bam,../01.mapping/S2.sort.rmdup.bam -S S1.splitters.bam,S2.splitters.bam -D S1.discordants.bam,S2.discordants.bam -o all.sv.lumpy.vcf

# step4_ind_genotype
# 上面lumpy只是在群体水平鉴定是否存在SV,但是没有鉴定具体样本SV位置
conda activate py2 # svtyper需要在python2.7的环境
for i in `seq 2`;do
        vcftools --vcf all.sv.lumpy.vcf --indv S${i} --recode --recode-INFO-all --out S${i} && svtyper -i S${i}.recode.vcf -B ../01.mapping/S${i}.sort.rmdup.bam -o S${i}.genotype.vcf
done

# step5_combine_vcf
# 将上一步的单个样品的genotype.vcf文件合并为一个文件列表vcf.list
ls S*.genotype.vcf > vcf.list
# 将genotype.vcf文件合并
svtools vcfpaste -f vcf.list > all.genotype.vcf

```
### CNV calling
- 软件：CNVnator
- 原理：read depth
- 运行过程：
    - 准备基因组文件目录
    - 运行cnvnator
```
$ cat CNV_pipeline.sh
#!/bin/bash
# cnvnator检测CNV

export ROOTSYS=/pub/software/root
export PATH=/pub/software/root/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib

cnvnator="/pub/software/CNVnator/cnvnator"
getvcf="/pub/software/CNVnator/cnvnator2VCF.pl"

ln -s ../data/genome.fasta genome.fa

# cnvnator需要基因组文件，且每条序列存储成一个文件且需要id.fasta这样命名
# 就是fasta文件中>跟的名称加.fasta
seqkit split -i genome.fa -O split  # 按id分割
ls ./split/genome.id_*.fa | sed 's/\.\/split\/genome\.id_//g' | awk '{print "mv ./split/genome.id_"$1" ./split/"$1}' > mv_name.sh

bash mv_name.sh

for i in `seq 2`;do
        bamfile="../01.mapping/S${i}.sort.rmdup.bam"

        # 将比对文件输出成root文件
        $cnvnator -genome genome.fa -root S${i}.root -tree $bamfile

        # 生成深度分布，每个窗口的
        $cnvnator -genome genome.fa -root S${i}.root -his 500 -d split

        # 统计计算，临近的窗口的深度差异
        $cnvnator -root S${i}.root -stat 500

        # 检查bin_size是否合适：作者推荐per bin中两个数相除在4-5之间
        $cnvnator -root S${i}.root -eval 500 > S${i}.eval.ratio

        # RD信号分割
        $cnvnator -root S${i}.root -partition 500

        # CNV检测
        $cnvnator -root S${i}.root -call 500 > S${i}.cnv

        # 格式转换
        $getvcf S${i}.cnv > S${i}.vcf
done
```
## step2_vairant_annotation
- 变异注释
    - 目的：查看变异相对于基因的位置；查看对基因结构，氨基酸编码的影响
- 输入文件
    - vcf格式变异结果
    - fasta
    - gtf，gff
- 软件：
    - annovar：
        - 构建数据库：
        - 注释：
            - 生成表格文件
    - snpeff：
### SNP_InDel
#### annovar注释SNP和INDEL
```
ar$ cat annovar.sh
#!/bin/bash

annovar="/pub/software/annovar"
genome_gtf="../data/genome.gtf"
genome_fasta="../data/genome.fasta"
snp_vcf="../data/all.filtered.snp.vcf"

# gtf => refGene(基因组浏览器常用)
gtfToGenePred -genePredExt $genome_gtf genome_refGene.txt
# gff3ToGenePred genome.gff genome_refGene.txt

# annovar基于refGene对序列进行提取序列，相当于把每个转录本提取出来
perl $annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile $genome_fasta genome_refGene.txt --out genome_refGeneMrna.fa

# 将vcf文件转换成表格形式
perl $annovar/convert2annovar.pl -format vcf4 -allsample -withfreq $snp_vcf > annovar.input

# 进行变异注释，如果需要所有信息添加 -separate
perl $annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver genome -dbtype refGene -outfile all.anno -exonsort annovar.input ./

```
### SV_CNV
## 群体结构分析
- 群体结构分析标配（基于SNP数据）
    - 进化树
    - PCA
    - STURCTURE

> 进化树(phylogenetic tree)：系统发育树，系统发育，是物种间，基因间，群体间乃至个体间谱系关系的一种表现形式：










