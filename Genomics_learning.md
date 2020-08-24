[TOC]
# 准备Loading ====>
## 体系图
**基因组学的分析**
- 基因组 survey
> 做一个物种基因组之前，初步估计一下这个物种的复杂情况（低深度的二代测序就可以搞定），基因组大小，基因组杂合度，基因组重复序列的含量，根据基因组的情况制定一个合适的测序策略（二代三代，光学图谱，HiC）
- 基因组的拼接
> 拿到各种测序数据之后，拼接基因组 => 物种的基因组序列
- 基因组结构注释
> 每个染色体什么位置是基因，什么位置是重复序列，然后基因上从哪到哪是外显子，从哪到哪是内含子
- 基因功能注释
> 每个基因是干嘛的，在通路上的作用
- 基因组评价
> 对基因组的拼接，注释做一个总体评价

**比较基因学分析（比较物种的基因组）**
- 物种层次
> 进化树构建与物种分化时间估计
- 基因组层次（基因组序列的改变）
> 共线性分析：物种分化之后基因组会发生变化，不同物种之间基因序列的变化；全基因组复制事件（WGD）研究：基因组加倍的事件，二倍体演变成四倍体，四倍体逐渐演变丢失又变成二倍体
- 基因家族层次
> 基因家族的鉴定：基因组上的所有基因进行家族化的划分（所有动植物源自于共同的祖先）；基因家族的收缩与扩张：哪些基因家族是某一个物种特有的=> 有和无的差异；哪些基因家族在某个物种多在另一个物种少=> 多和少的差异
- 基因层次：
> 正选择分析：人和大猩猩有一个基因是直系同源基因，但是这个基因序列在人和大猩猩中不一样，就可以去研究是否这个基因受到了正选择，是什么导致的，环境？驯化？
## 关于测序
通量：一开机能产生的数据多少，意味着测序的价格(1024是计算机存储问题，生信10进制)
- 1k = 1000bp
- 1M = 1000 000bp
- 人的基因组数据大概3G，测了90G数据，相当于把基因组的每一个位置测了30遍（深度）

主流测序平台：
- **illumina**：边合成边测序，读长150bp，通量6Tb
    - dNTP的荧光标记 => 区分不同碱基
    - 桥式PCR形成cluster => 加强荧光信号
    - 末端终止法，暂停合成 => 识别荧光信号
| 常见 | 含义 |
| --- | --- | 
| Read: | 测序得到的片段，read length：能测的序列长度 |
| inset size | 待测序片段的长度 |
| 测序深度depth | 每个位置测了多少次，随机性，深度=测序量/基因组大小（DNA片段打断 => 文库）|
| single/paired end | 单端/双端测序 |
| coverage | 测序覆盖率 | 
- PacBio：边合成边测序，读长14k，通量20Gb
    - 荧光标记，碱基在聚合酶的时间越长荧光越强，illumina是在一个很空旷的空间识别荧光而PacBio是在一个很狭窄的小孔空间识别，所以不需要荧光放大
    - PacBio的接头为哑铃型（illumina的为Y型），优点：双链打开后是一个环形，边合成边测序之后扔掉一个环，重新再测序实现对一条链的反复测序提高准确率；无GC偏好（GC含量高PCR效率低）
- Nanopore：电流信号，读长>2M，通量3~6Tb，牛逼就完事
## 软件安装
### 基因组学重要的软件
| 名称 | 功能 | 干啥 |
| --- | --- | --- | 
| RepeatMasker | 重复序列屏蔽 | 已知哪些序列是重复序列，在基因组中找到屏蔽掉，换成N | 
| RepeatModeler | 重复序列鉴定 | 找一下基因组上都有哪些重复序列 | 
| Braker | 基因参数模型训练 | 太复杂了总之，每一个物种的基因偏差是不一样的 |
| Maker | 基因结构注释 |
| MCScanX | 共线性分析 | 
| GCBI | 共线性分析，一个新的python开发的 |
### RepeatMasker安装：重复序列搜索 
```
$ conda install RMBlast TRF #RepeatMasker依赖的两个包, RMBlast专门为repeatmasker开发的搜索引擎,用来取数据库的
#还需要重复序列的数据库：Repbase（各种物种的重复序列收集整理放在数据库中,因此找到的重复序列都是已知的,位置的可以参考RepeatModeler）
#重复序列：简单重复序列：卫星序列，转座子（LTR，RNA转座子）
#下载RepeatMasker.tar，解压后进入文件，将Repbase的压缩文件放入repeatmasker的library下
#查看README，得知查看INSTALL，得知运行：
$per ./configure # 首先是配置，一般Perl语言不需要配置，但是对于复杂的软件包需要，我感觉是因为要选择一些搜索引擎导入数据库的原因，配置完之后就会得到一个perl的可执行文件
#然后输入trf的位置，选择RMBlast搜索引擎，5.done        
$ perl RepeatMasker #配置完成后得到perl的可执行文件
```
### RepeatModeler安装：重复序列鉴定
```
$conda install recon RepeatScout Ltr_retriever MAFFT CD-HIT #同样先用conda装能装上的，RepeatModeler依赖的包很多

# 下载GenomeTools，解压（从github上或者官网上下载，安装比较麻烦，参数很多，下面的别人的经验）
$make 64bit=yes thread=yes cairo=no -j4 #编译的过程，将c语言编译成2进制
$make 64bit=yes thread=yes cairo=no -j4 install prefix=~/gt #拷贝到制定位置，不加install则默认拷贝到当前

#NINJA 安装，这个也在github下载最新的安装包，解压，找readme，发现用$make all命令编译
$make all

#终于可以安装RepeatModeler了，下载好压缩包解压，看一眼readme
$ perl ./configure #开始配置，应该是对程序进行设置，配置它的搜索引擎啊，数据库啊之类的，总之conda装的软件如果程序识别到了就不许改，手动装的程序就需要手动填写程序位置
#1.选择perl的位置，就用conda里的；
#2.选择RepeatMasker的路径（RepeatModeler依赖RepeatMasker），
#神坑，不能用conda里的（为啥子conda里面有呢，可能是装其他软件有的依赖，所以装了，但是不用这个，用我们配置好的），进入手动安装的RepeatMasker文件夹，找到可执行文件pwd获得路径 ；
#3.recon，conda装的跳过；
#4.rscout（repeatscout），conda装的跳过；
#5.trf，用which找到；6.选择RMBlast；7.Done；
#8.新版的repeatmodeler支持LTR的检测，选择配置，输入genometools的路径，一定要找到可执行文件，bin/下；
#9.Ltr_retriever，MAFFT，conda安装的跳过；10.NINJA，手动配置；11.CD-HIT，conda装的跳过；
 $ perl ./RepeatModeler #看一下是否安装成功，成功渡劫

#编译小知识：一个c语言编写的程序可以直接用gcc进行，但是一个软件包里面有很多个程序需要
#编译，为了提高编译效率加快编译速度，因为程序间有的软件有依赖关系所以需要按顺序编译，所以
#程序里边一般有一个Makefile文件，里面记录哪些软件依赖哪些软件，输入make就会按照Makefile
#里的顺序编译，所以make一定要在有Makefile的文件夹中运行，编译过程有些参数也需要调整
#具体看需求，编译完成后有可能报错，就需要清理掉之前的编译结果重新编译，因为编译后的
#程序会自动跳过，输入make clean删除刚刚的编译，就可以重新编译
```
### MCScanX的安装：共线性分析
`msa.h`,`dissect_multiple_alignment.h`,`detect_collinear_tandem_arrays.h`中的第一行添加`#include <unistd.h>`
```
# MCScanX上古软件，作者用32位写的所以需要改源代码支持64位，查看readme
```
### BRAKER安装：基因参数模型训练
超级难安装，BRAKER依赖的包尤其的多，而且最坑的是大部分用conda安装的都不能用，这个是别人的经验，我觉得如果是这样的话最好所有软件都手动安装，下面的这个大哥依赖的包
- Braker
    - augustus（基因注释）
        - samtools
        - bamtools
        - htslib
        - bcftools
        - tabix 
    - GeneMark-ET（基因预测）
    - GenomeThreader 
```
# BRAKER依赖的包
- AUGUSTUS 3.3.3
- GeneMark-ES/ET/EP 4.48_3.60_lic
- BAMTOOLS 2.5.1
- SAMTOOLS 1.7-4-g93586ed
- ProtHint 2.4.0
- GenomeThreader 1.7.0
- Spaln 2.3.3d 
-  (Exonerate 2.2.0 )
- NCBI BLAST+ 2.2.31+ 
- DIAMOND 0.9.24
- cdbfasta 0.99
- cdbyank 0.981

#依赖的perl包
File::Spec::Functions
Hash::Merge
List::Util
Logger::Simple
Module::Load::Conditional
Parallel::ForkManager
POSIX
Scalar::Util::Numeric
YAML
Math::Utils
# cpanm是perl语言中模块的管理工具
$cpanm File::File::Spec::Functions Hash::Merge \
List::Util Logger::Simple \
Module::Load::Conditional \
Parallel::ForkManager \
POSIX Scalar::Util::Numeric YAML Math::Utils

$conda install -y spaln exonerate #先用conda装，有一些conda装的软件不能用，就很奇怪，-y直接装不要问

# 安装Bamtools
# 下载解压后：查看readme 
# CMakeLists.txt 一般用cmake，cmake比make高级一点点，
$ mkdir ~/jimmy/software/bamtools #创建一个文件夹，软件安装在这里
$ cd /path/to/bamtools
$ mkdir build #专门用于编译软件
$ /cd build
$ cmake -DCMAKE_INSTALL_PREFIX=~/jimmy/software/bamtools .. #其实是在配置，并指定位置
$ make # 编译
$ make install # 安装，本质是将当前的拷贝到指定的位置

# 安装augustus：基因注释软件，从头预测，不依赖任何东西直接判断哪些是基因哪些是内含子外显子
# 需要修改augutus源代码因为可能会包bamtools的错，需要更改它所依赖的bamtools的路径
$ vim auxprogs/bam2hints/Makefile
11 INCLUDES = /home/u2052/jimmy/software/bamtools/includes/bamtools
12 LIBS = /home/u2052/jimmy/software/bamtools/lib/libbamtools.a -lz
$ vim auprogs/filterBam/src/Makefile
11 BAMTOOLS = /home/software/src/bamtools/include/bamtools
12 INCLUDES = -I$(BAMTOOLS) -Iheaders -I./bamtools
13 LiBS = /home/u2052/jimmy/software/bamtools/lib/libbamtools.a -lz 
# augustus依赖的包如下：
# (samtools htslib bcftools)这三位神仙要放在同一个文件夹下以软件名为名  tabix
# 安装这四个去逛网下或者github优先考虑release
# 注意tar.bz2的解压命令 tar -jxvf xxx.tar.bz2
# 修改bam2wig的Makefile指向TOOLDIR指向以上安装路径
$ vim auxprogs/bam2wig/Makefile
9 ifndef TOOLDIR
10     TOOLDIR=/home/u2052/software
# => 完成上述，在augustus目录下make
$ cd augustus-3.3.2
make

#将augustus添加到环境变量
export PATH=/home/u2052/jimmy/software/augustus-3.3.2/bin:$PATH
export PATH=//home/u2052/jimmy/software/augustus-3.3.2/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=/home/u2052/jimmy/software/augustus-3.3.2/config

# 下载GeneMark-ET，官网需要填写信息，同时记得下载证书，证书解压后放在home目录下改名为 .xxx
# 为perl语言编写的直接可以用，但是注意当前有很多perl，需要换成我们当前的perl因为，在此之前我们对perl加载了很多模块，用下面命令批量换
ls *.pl | xargs -i sed -i 's/#!\/usr\/bin\/perl/#!\/usr\/bin\/env perl/' {}

# 装大哥大BRAKER
# 下载后加压完事，添加一个环境变量
$ export PATH=/home/u2052/jimmy/software/BRAKER-2.1.2/scripts:$PATH
```
### Maker：基因注释
基因注释方法：从头预测用est用蛋白去做RNA-Seq去做=> 整理证据数据 =>格式不够标准上传到数据库啥的会被质疑
> orthofinder：-py2.7生信基因功能分析工具
## 任务：
### 任务一：数据过滤和质控
```
fastp -i Col-0-0h-rep1_1.fq.gz -I Col-0-0h-rep1_2.fq.gz -o ../clean_data/Col-0-0h-rep1_1.fq.gz -O Col-0-0h-rep1_2.fq.gz -h Col-0-0h-rep1.html -j Col-0-0h-rep1.json 1> ../clean_data/Col-0-0h-rep1.log 2>&1 &
```
### 任务二：基因家族分析
> background：芝麻在一些油脂方面的合成比其他作物丰富，其可能原因基因家族的数目层次比较多（基因扩展），也有可能是表达上的，研究FAD这个基因家族，看一看这个家族在芝麻中多少个拷贝，在大豆中多少个拷贝，拟南芥研究的比较透彻有3个，于是就用拟南芥的这三条基因作为query序列，拿大豆和芝麻的蛋白序列作为database

 - 输入数据：query
    - 拟南芥的3条FAD4-like，蛋白序列（优先选择蛋白序列，因为蛋白序列更保守一些））
    - 芝麻的蛋白集合
    - 大豆的蛋白集合
    - 拟南芥的蛋白集合
- 分析软件：blast
- 分析步骤：
```
#构建blast数据库
$ cat atha.fasta gmax.fasta sind.fasta > all.fasta
$ makeblastdb -dbtype prot -in all.fasta -out ../blastdb/prot

#blast比对
$ blastp -query data/atha_FAD4.fa -db blastdb/prot -out result.txt #得到一个详细的报告适用于几个样，默认Pairwise
$ blastp -query data/atha_FAD4.fa -db blastdb/prot -out result.txt -outfmt 7 -evalue 1e-10 #推荐这个格式，蛋白质比对evalue一般设置1e-5到1e-30之间
$ less -S result.txt
# BLASTP 2.9.0+
# Query: atha|AT1G62190.1
# Database: blastdb/prot
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, eval
# 10 hits found
atha|AT1G62190.1        atha|AT1G62190.1        100.000 295     0       0   
atha|AT1G62190.1        atha|AT4G27030.1        68.770  317     70      2       
atha|AT1G62190.1        atha|AT2G22890.1        70.775  284     76      4       
atha|AT1G62190.1        gmax|KRG98807   54.861  288     121     
atha|AT1G62190.1        gmax|KRH46136   55.160  281     106     
atha|AT1G62190.1        sind|SIN_1014454        50.336  298     120     4     

# 对result.txt进行过滤
$ grep '^[^#]' result.txt | awk '$3 >= 50 && $11 <= 1e-30{print $2}' | sort | uniq > gene_list.txt #筛选
$ grep -v '#' result.txt 
$ cat gene_list
atha|AT1G62190.1
atha|AT2G22890.1
atha|AT4G27030.1
gmax|KRG98807
gmax|KRH46136
sind|SIN_1005002 
sind|SIN_1005003
sind|SIN_1005004
sind|SIN_1014454

# 根据gene_id提取序列 ：seqtk， seqkit
$ seqtk subseq data/all.fasta gene_list.txt > FAD4.fasta
$ head FAD4.fasta
>atha|AT1G62190.1
MAVSFQTKNPLRPITNIPRSYGPTRVRVTCSVTTTNPQLNHENLVVEKRLVNPPLSKNNDPTLQST
WTHRLWVAAGSTTIFASFAKSIIGGFGSHLWLQPALACYAGYVFADLGSGVYHWAIDNYGGASTPI
VGAQLEASQGHHKYPWTITKRQFANNSYTIARAITFIVLPLNLAINNPLFHSFVSTFAFCILLSQQFH
AWAHGTKSKLPPLVMALQDMGLLVSRKDHPGHHQAPYNSNYCVVSGAWNKVLDESNLFKALEM
ALFFQFGVRPNSWNEPNSDWTEETETNFFTKI
>atha|AT2G22890.1
MATSLQTKYTLNPITNNIPRSHRPSFLRVTSTTNSQPNHEMKLVVEQRLVNPPLSNDPTLQSTWTH
RLWVAAGCTTVFVSFSKSIIGAFGSHLWLEPSLAGFAGYILADLGSGVYHWATDNYGDESTPLVGIH
IEDSQDHHKCPWTITKRQFANNLHFMARGTTLIVLPLDLAFDDHVVHGFVSMFAFCVLFCQLFHAW
AHGTKSKLPPLVVGLQDIGLLVSRIHHMNHHRAPYNNNYCVVSGVWNKVLDESNVFKAMEMVLYIQ
LGVRPRSWTEPNYE
```
# 基因组学：
## 主要内容：
- 基因组拼接
- 基因组注释
    - 基因注释
    - 重复序列注释 
- 比较基因组学分析
    - 进化树构建 
    - 基因家族聚类
    - 全基因组复制
    - 共线性分析
## => STEP1_基因组组装（De novo assembly）：
> 对一个序列未知的物种进行测序，从头构建其参考基因组的过程

### 数据质控（survey之前）：
- 过滤接头，N含量，低质量
- 检测GC含量及分布
- 去除duplication（可选）
- NT库比对（如果GC含量分布很奇怪，可以看一看是什么物种导致的污染）
- FastQC
### 基因组组装和注释理论：
参考润泽的简书：[基因组组装和注释理论]<https://www.jianshu.com/p/031de137bf38>
### 01.基因组survey：
Survey用少量的数据评估基因组基本信息的有效手段，对于没有参考基因组的物种，对基因组信息有一个明确的概念，对后续的测序及组装方案是很有必要的（**一般测序量为预估基因组大小的30-50X，二代测序**）
- 通过少量的测序数据的survey可以知道一下信息：
    - 基因组大小：决定测序量和组装深度
    - 重复序列比例
    - 杂合度：对于二倍体，两组染色体会有差异
    - GC含量：过高或过低的GC含量会导致测序偏向性，影响基因组分析结果
#### jellyfish的使用：kmer的频数表
```
$ ~/jimmy/workspace/01.Assembly/01.kmer_analysis$ cat step1.jellyfish.sh
#!/bin/bash 
# jellyfish统计kmer频数表
jellyfish=/home/u2052/miniconda3/bin/jellyfish

#jellyfish设置kmer值
pre=Kmer_19


#jellyfish必须是解压好的文件
ls /home/u2052/jimmy/workspace/data/fastq/ERR1938683*.fastq.gz | awk '{print "gzip -dc "$0}' > generate.file

#step1_count基于fastq中的reads做kmer统计，kmer的碱基序列，在序列中出现的频数，生成一个二进制文件
$jellyfish count -t 4 -C -m 19 -s 1G -g generate.file -G 2 -o $pre
# -t :线程数
# -C :reads的正负链都考虑
# -m :kmer值
# -s :初始的内存大小
# -g :包含解压命令的文本
# -G :解压文件的线程数
# -o :输出

#step2_histo继续count中的碱基序列和频次做统计得到kmer频数表，|kmer的深度：kmer出现的次数|种类数有多少|
$jellyfish histo -v -o $pre.histo $pre -t 4 -h 10000
# -v：生成日志文件
# -o：输出
# ：输入
# -h：统计深度的阈值

#step3_stats：kmer的总个数，kmer的总类数，kmer的uniq的次数（kmer出现次数为一的个数）
$jellyfish stats $pre -o $pre.stat
```
#### gce ：基于kmer频数表算基因组大小
```
$ cat gce.sh
#!/bin/bash
/pub/software/gce-1.0.2/gce -f Kmer_19.histo -c 60 -H 1 -g 876075288 -M 10000 > gce.table 2>gce.log
# -f：输入kmer频数表
# -c：峰值的位置
# -H：杂合设置1，非杂合设置0
# -g：总的kmer个数
# -M：考虑到的最大深度
# 最终需要的文件再gce.log里面迷幻操作
```
#### GenomeScope ：基于kmer频数表算基因组大小
（基于R，比实际基因组大小小一点）
```
$ Rscript /pub/software/genomescope/genomescope.R Kmer_19.histo 19 150 ./ 100000
# jellyfish生成的kmer频数表 
# kmer数
# reads长度
# 分析的kmer最大有效深度K
# 两个plot图和一个summary文件
```
> 重点：
杂合情况：主峰深度1/2出现峰值
重复情况：主峰2倍左右出现峰值
出现多个主峰：有可能存在污染
### 02.基因组组装：
reads => contig => scaffold => chromosome  
- => contig组装：（连续的没有N的序列，二代基于kmer，三代overlap）
    - 算法：
        - OLC: Over lap Layout Consensus（三代）
            - **Canu**(准确度很高，但是很慢，一般保留杂合区域） 
            - **WTDBG2**（速度极快，但对重复度高的基因组准确度不够，因为 会忽略掉杂合区域）
            - **Flye** 
            - Celera Assembler
            - Falcon 
            - Miniasm
        - DBG: De Brui jn Graph（二代）Kmer
            - **SOAPdenovo2**（拼接速度快，直接conda）
            - **SPAdes** （适合小基因组尤其是200M以下，消耗资源很多）
            - AllPath-LG
            - Velvet
- => scaffold组装：基于contig的前后关系再次连接，中间位置的序列插入N，是contig的长度更长，如果contig指标不好可以构建，如果好的话可以跳过。三代一般不构建 
    - **mate pair**
    -  **bionano**
    - 10x genomics
- => polish
    - pilon.jar
    - 二代测序数据
    - 30x-100x
- => chromosome：
    - 遗传图谱挂在染色体（**ALLMAPs**）
    - Hi-C连接染色体i（**ALLHiC**，**3d-DNA**）
> Hi-C：同一条染色体上的位置更容易被蛋白交联在一起（空间距离个更近），测被蛋白交联的序列，跟据这个序列的数目，将contig分组到各个染色体

> 组装结果评价：
N50：？
NG50：？
#### soapdenovo2：2代组装
```
$ cat config.txt
max_rd_len=150

[LIB] #如果有多组数据写多个[LIB]
avg_ins=250 #文库片段为250
reverse_seq=0 #非环化文库为0，环化文库为1
asm_flags=3 #即参与contif的拼接，又参与scaffold的构建
rd_len_cutoff=150 #reads的长度
rank=1 #
pair_num_cutoff=3 #有多少个reads支持，那么可信，一般为3，5，7
map_len=32 #
q1=/home/u2052/jimmy/workspace/data/fastq/ERR1938683_1.fastq.gz
q2=/home/u2052/jimmy/workspace/data/fastq/ERR1938683_2.fastq.gz


$ cat run_soapdenovo.sh
#!/bin/bash
SOAPdenovo=/home/u2052/miniconda3/envs/soapdenovo2/bin/SOAPdenovo-63mer

SOAPdenovo pregraph -s ./config.txt -o K41 -K 41 -p 8 -d 1 2>./pregraph.log
# -o：指定输出文件的前缀
# -K：拼接的kmer大小
# -d：kmer出现频率为1的删掉
SOAPdenovor contig -g K41 -p 8 2>./contig.log
SOAPdenovo map -s ./config.txt -g K41 -k 35 -p 8 2>./map.log t
# -k：kmer选择35：reads更多的比对到contig上，kmer尽可能选择奇数
SOAPdenovo scaff -g K41 -p 8 -F 2>./scaff.log
# -F：补洞

$ bash run_soapdenovo.sh
$ ll
K41.scafSeq #拼接的结果
K41.scafStatistics # 统计
```
#### sapdes：2代组装
（之前是对小基因组进行拼接200M以下，非常消耗资源）
```
python /pub/anaconda3/bin/spades.py --careful -t 6 -k 51,127 --pe1-1 /home/u2052/jimmy/workspace/data/fastq/ERR1938683_1.fastq.gz --pe1-2 /home/u2052/jimmy/workspace/data/fastq/ERR1938683_2.fastq.gz -o spadeD
```
#### canu：三代拼接
```
path/to/canu -pacbio /home/u2052/jimmy/workspace/data/fastq/DRR076709_RSII.fastq.gz -p test -d ./output genomesSize=14m useGrid=false
```
#### miniasm：三代拼接
（不建议用，如果数据质量不好不出结果）
#### wtbdg2：三代拼接，速度超快
```
# 从github上下载之后解压，直接生成可执行文件
export PATH=/home/u2052/jimmy/software/wtdbg-2.5_x64_linux:$PATH

$ cat run_wtdbg2.sh
#!/bin/bash
# 组装
wtdbg2 -t 6 -x rs -g 14M -L 500 -l 500 -e 2 -i /home/u2052/jimmy/workspace/data/fastq/DRR076709_RSII.fastq.gz -o test
# -x：数据类型
# -g：基因组大小
# -L：使用的最小reads长度
# -l：reads之间overlap的最小长度
# -e：结果中支持度小于2的不要，默认是3
# 拼接的结果为test.ctg.lay.gz

# 得到一直性序列fasta
## wtdbg-cns 2分钟
wtdbg-cns -t 6 -i test.ctg.lay.gz -f -o test.wtdbg-cns.fa
## wtpoa-cns12分钟
wtpoa-cns -t 6 -i test.ctg.lay.gz -f -o test.wtpoa-cns.fa
# -f：若文件存在，终止，-f强制覆盖原文件
# 转化后的文件test.wtdbg-cns.fa
```
#### pilon：三代组装的拼接矫正
（二代数据矫正三代拼接结果）
> 二代单碱基错误率叫少可以跳过矫正这一步，也可以三代对三代进行矫正，用的是recon，但没必要，因为三代矫正效果不明显

```
$ cat run_pilon.sh
#!/bin/bash 
draft_genome=/home/u2052/jimmy/workspace/data/genome/genome.fasta #三代组装结果
fq1=/home/u2052/jimmy/workspace/data/fastq/ERR1938683_1.fastq.gz
fq2=/home/u2052/jimmy/workspace/data/fastq/ERR1938683_2.fastq.gz
samtools=/home/u2052/miniconda3/bin/samtools

ln -s $draft_genome ./genome.fa
bwa index genome.fa
bwa mem -t 6 genome.fa $fq1 $fq2 | $samtools sort -o reads.bam
$ samtools index reads.bam

java --Xmx10G -jar /pub/software/pilon-1.23/pilon-1.23.jar \
-- genome $draft_genome \
-- frags reads.bam \
--changes --vcf --diploid --threads 6 \
--outdir ./pilon_out --output genome_pilon
```
> assembly-stats：对基因组进行简单的统计，直接用conda安装，可以得到n50 - n90的指标
### 03.Hi-C进行染色体挂载（使用ALLHi-C）
需要指定染色体条数，allhic会按照染色体的条数进行聚类，把contig聚成染色体，然后排序
- 输入：
    - 染色体拼接结果（contig，scaffold，非染色体水平的拼接结果）
    - Hi-C测序数据（PE150）
    - 染色体条数
    - 酶切类型
- Hi-C比对到contig：
    - bwa index
    - bwa aln
    - bwa sample
    - 保留酶切位点附近的且成对的比对结果：`PreprocessSAMs.pl`
    - 对比对的质量进行过滤：`filterBAM_forHiC.pl`
- allHic连接：
    - contig聚类成cluster => .txt：
        - 需要知道染色体条数，酶切类型
        - 比对的bam文件
        - 组装的基因组
    - .txt => CLM文件：酶切位点个数等信息
- 排序和方向优化 => .tour：
    - sample.clean.clm
    - xxx.5g1 g2
- `ALLHiC_build`生成染色体：输入基因组文件即可，会自动识别.tour => **groups.agp** 和 **groups.asm.fasta**
- heatmap：`ALLHiC_plot` + xxx.bam + 染色体长度文件
```
# bwa有两种比对模式，mem和aln（适合短reads，比对要求更严格）
$ cat run_allhic.sh
#/bin/bash
## 加载环境
export PATH=/pub/software/ALLHiC/scripts/:/pub/software/ALLHiC/bin/:$PATH

k=5
genome=./Ler-1.allpaths_lg.final.assembly.fasta
fq1=/home/wangying/RD/01.genome/Genome_wy/workspace/01.Assembly/00.data/NNJHiC2_1.fastq.gz
fq2=/home/wangying/RD/01.genome/Genome_wy/workspace/01.Assembly/00.data/NNJHiC2_2.fastq.gz

## 构建index
ln -s  $genome ./draft.asm.fasta
bwa index  draft.asm.fasta  
samtools faidx draft.asm.fasta  

## 比对
bwa aln -t 12 draft.asm.fasta $fq1 > sample_R1.sai  
bwa aln -t 12 draft.asm.fasta $fq2 > sample_R2.sai  
bwa sampe draft.asm.fasta sample_R1.sai sample_R2.sai $fq1 $fq2 > sample.bwa_aln.sam  # 合并

## 对比对结果进行过滤  HINDIII: AAGCTT , MBOI: GATC（Hi-C的文库构建有酶切，位点是确定的，利用这个对结果进行过滤）
PreprocessSAMs.pl  sample.bwa_aln.sam draft.asm.fasta  HINDIII # 过滤酶切位点和成对的
filterBAM_forHiC.pl  sample.bwa_aln.REduced.paired_only.bam sample.clean.sam # 过滤比对到多处的reads
samtools view -H sample.bwa_aln.bam > header
cat header sample.clean.sam | samtools view -bS - > sample.clean.bam #bam做为hic主程序的输入
# samtools view -bt  draft.asm.fasta.fai  sample.clean.sam > sample.clean.bam（这个版本的allhic有点问题，无法识别header，这一行本来是作者的代码，改用上面的步骤）

## 基于比对结果将contig划分成group
 ALLHiC_partition -b sample.clean.bam -r draft.asm.fasta -e AAGCTT -k $k

## Extract CLM file and counts of restriction sites 统计每个contig上的酶切位点个数，比对到contig的条数，等等数据，少于要求的过滤掉，结果用于排序，
allhic extract sample.clean.bam draft.asm.fasta --RE AAGCTT

## 对group内的contig排序和方向优化
for K in `seq 1 $k`
do 
    allhic optimize sample.clean.counts_AAGCTT.*g${K}.txt sample.clean.clm
done

## 得到最终结果
ALLHiC_build draft.asm.fasta

## 绘制heatmap 图

seqkit   fx2tab -l -i  -n  groups.asm.fasta > len.txt
grep 'sample.clean.counts_AAGCTT' len.txt > chrn.list
ALLHiC_plot sample.clean.bam groups.agp chrn.list 500k pdf
```
### 04.GC-depth
将组装好的基因组划分成等大的窗口，统计每个窗口平均的GC含量和平均的测序深度 => 判断是否有污染
如果发现有污染可以比对一下 NT核酸数据库（NR：NCBI非冗余蛋白数据库），或者16s数据库
```
$ cat step1.baw.sh
#!/bin/bash
# 测序的reads比对到组装好的基因组
ln -s ../00.data/genome.fasta
bwa index  genome.fasta
bwa mem -t 4  genome.fasta ../00.data/ERR1938683_1.fastq.gz ../00.data/ERR1938683_2.fastq.gz | /pub/software/samtools/samtools sort - -o aln_sort.bam

$ cat step2.gc_depth.sh
#!/bin/bash
genome=genome.fasta
bam=aln_sort.bam
prefix=gcdep
window=500
step=250

# 计算序列长度
seqtk comp  $genome  | awk '{print $1"\t"$2}' > $prefix.len # seqtk comp：统计每个序列的碱基个数碱基含量等很实用

# 划分窗口 生成bed文件
bedtools  makewindows -w $window -s $step -g $prefix.len > $prefix.window.bed 
# -s：滑动窗口 ，1-500 250-750 ... （从0开始计数，左闭右开[) ，我感觉是左开右闭）

# 按窗口提取序列并计算gc含量  
seqtk subseq $genome  $prefix.window.bed  > $prefix.window.fasta
seqtk comp  $prefix.window.fasta |awk '{print $1 "\t" ($4+$5)/($3+$4+$5+$6) } ' > $prefix.window.gc

# 按窗口计算平均深度
bedtools  coverage -b aln_sort.bam -a gcdep.window.bed -mean  | awk '{print $  1":"$2+1"-"$3"\t"$4}' >  $prefix.window.depth
# coverage：统计深度信息
# -b : 比对好的bam文件
# -a：窗口的bed文件
# -mean：以窗口为单位计算平均深度

# 绘图
Rscript  run_gcdep.R $prefix.window.gc $prefix.window.depth $prefix.pdf 0 0.8 0 500

$ cat run_gcdep.R
#!/bin/rscript
args<-commandArgs(TRUE)

if(length(args) < 3 ) {
  stop("USAGE: run_gcdep.R  gc.table  depth.dable  out.pdf  min_gc  max_gc  min_dep max_dep \n")
} else{
  file1 = args[1]
  file2 = args[2]
  outpdf = args[3]
  mingc = as.numeric(args[4])
  maxgc = as.numeric(args[5])
  mindp = as.numeric(args[6])
  maxdp = as.numeric(args[7])
}

args[6]
args[7]

library(tidyverse)
gc <- read.table(file1, header = F , row.names = 1)
dp <- read.table(file2, header = F , row.names = 1 )

gc1 <- rownames_to_column(gc,var = "id")
dp1 <- rownames_to_column(dp,var = "id")

gcdep <- inner_join(gc1, dp1, by = "id")

colnames(gcdep) <- c("id", "gc", "dep")

p <- ggplot( data = gcdep )  +
  geom_point(aes(x=gc,y=dep), colour = "red", alpha = 0.5 ,cex = 0.7) +
  xlab("GC content")  +
  ylab("Depth") +
  xlim(mingc, maxgc) +
  ylim(mindp, maxdp) +
```
## => STEP2_基因组注释：
（所有软件皆可用自己安装的，部分代码没有更改路径）
- 重复序列注
    - => 简单重复序列（串联重复序列，卫星序列）（Tandem repeats）：亲子鉴定，遗传标记（主要分布在着丝粒区域）
    - => 转座子（散在重复序列）：  
        - DNA转座子:ClassII
            - （subclass）TIR
            - Crypton
        - RNA转座子:Class I
            - **LTR**:long terminal repeat:
                - 计算LTR的插入时间：当一个LTR转座子刚刚插入到基因组上的时候，左右两个LTR是一摸一样的，随着时间的演化，这两个LTR分别突变分别演化最后越来越不像，分子距离S（描述他们之间的差异）越来越远，S/进化的速度t => LTR大概的插入时间。算出多个LTR的插入时间，得到一个分布，一般期望位置就是LTR的插入时间。注意：突变的速度不是恒定的受到自然和人为的选择，因次LTR的选择位点很关键，选择同义替换位点。
            - DIRS
            - PLE
            - **LINE**
            - **SINE**

>（比较基因组学）Whole genome duplication（物种在演化过程中基因组的加倍，长期以后逐渐又丢失）
>（基因家族）基因的多拷贝=> 大的基因家族

- 基因结构注释
- 蛋白功能注释
### 重复序列注释：
- 注释原理：
    - 已知重复序列：（大概找一找）
        - Dfam（只包含DNA转座子）
        - Repbase（大多数物种重复序列数据库，比对一下他）
    - 未知重复序列：
        - 简单重复序列：`trf`
        - 重复次数较高的序列：`RECON`,`RepeatScout`（都是基于kmer）
        - LTR（比例相当大占基因组的20%）：`LtrHarvest`,`Ltr_retriever`
#### 01.构建RMBLAST数据库：
> 提取某条染色体的序列：$ seqtk subseq genome.fasta sample.txt > subgenome.fasta 

> seqtk也支持bed文件 LG1    0    1000 (提取第一个到1000个)

> $ fasta_formatter -i subgenome.fasta -o subgenome_format.fasta -w 60 # 一行换成60行`FASTAX-Toolkit`

- 输入：基因组序列（genome.fasta）
- 输出：RMBLAST数据库（genome.nhr等）
- 命令：`BuildDatabase`
#### 02.生成物种自身的重复序列数据库：
- 输入：RMBLAST数据库(genome.nchr等)
- 输出：物种自生重复序列数据库(sp-families.fa)
- 命令：`RepeatModeler` => 非冗余数据库(包含了所有的trf,recon,repeatscout,ltrharvest,ltr_retriever)       
#### 03.鉴定和屏蔽重复序列
(具体到每个位置有哪些重复序列)
- 输入：
    - 基因组序列(genome.fasta)
    - 物种自身重复序列数据库(sp-families.fa)
- 输出：
    - 统计结果(genome.fasta.tbl)
    - 详细重复序列信息(genome.fasta.cat)
    - 表格格式重复序列信息(genome.fasta.out)
    - 屏蔽后的基因组(genome.fasta.masked)
- 命令：`RepeatMasker`
#### 04.画图：
- 输入：详细重复序列信息(genome.fasta.cat)
- 输出：RepeatLandscape(RepeatLandscape.html)
- 命令：`calcDivergenceFromAlign.pl``createRepeatLandscape.pl`
```
$ cat run_repeat.sh
#!/bin/bash
# 构建RMBLAST数据库
/home/u2052/jimmy/software/RepeatMasker/BuildDatabase -name sesame test_data/genome.fasta
# 输出=> html nhr nin nnd nni nog snq 等等

# 运行 RepeatModeler构建组装好的构建重复序列数据库
/home/u2052/jimmy/software/RepeatMasker/RepeatModeler -database sesame -pa 20 -LTRStruct
# -pa 20：线程数
# 输出：
# stk文件，种子文件
# families.fa :物种的重复序列数据库，即包含已知的又包含特有的

# 运行 RepeatMasker组装的基因组 + 这个基因组重复序列数据库 = 重复序列的具体位置
/home/u2052/jimmy/software/RepeatMasker -pa 20 -qq -lib sesame-families.fa test_data/genome.fasta >repeatmasker.log 2>&1
# -qq：超快
# -s：慢速
# -q：快
# 输出：输出到基因组所在的文件夹
# .cat ：超详细没人看，除非研究重复序列的科学家
# .out：有点像bed文件，一般也不看
# .tbl：类似于一个统计报告，一般要放文章
# .fasta.masked：屏蔽之后的

# 生成RepeatLandscape
gunzip genome.fasta.cat.gz
# 根据比对结果计算分子距离
perl /home/u2052/jimmy/software/RepeatMasker/util/calcDivergenceFromAlign.pl -s sesame.divsum test_data/genome.fasta.cat

perl /home/u2052/jimmy/software/RepeatMasker/util/createRepeatLandscape.pl -div sesame.divsum -g 18577337 > sesame.html
# -g：指定基因组大小
```
### 基因结构注释：(编码蛋白的区域)
> 屏蔽掉重复序列的方式：1.重复序列替换成小写的字母softmask；2.直接替换成大写的N，hardmask

- 原理：
    - **RNA层次证据**：尽可能取多种因素使其转录，挖掘更多信息（20个样本，20个时间点构成一个文库，深度100M以上，指的是序列的条数）
        - 数据源：EST序列，二代测序RNA拼接，三代测序全场转录本
        - 用法：blastn比对到基因组
    - **Protein层次证据**（一般来自于近缘物种，既近缘研究又清楚的）：
        - 数据源：数据库下载的蛋白序列
        - 用法：exonerate比对到基因组
    - ab-initio预测证据（从头预测）（机器学习算法：HMM）：
        - 参考模型：已有的物种参数（待注释物种的近缘物种），自训练参数模型（**高可信的数据源，比如上面两个层析的数据**）（`Braker`训练数据）
        - 用法augustus预测

=> 证据整合
- 软件：（**Maker**整合证据，GLEAN，EVM）
- 输出：（gff/gtf文件，cds序列，protein序列）（gff的exon序列=> cds序列 => protein序列)
```
# 添加一下环境变量
# augustus
export PATH=/home/u2052/jimmy/software/augustus-3.3.2/bin:$PATH
export PATH=/home/u2052/jimmy/software/augustus-3.3.2/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=/home/u2052/jimmy/software/augustus-3.3.2/config
# Braker
export PATH=/home/u2052/jimmy/software/BRAKER-2.1.2/scripts:$PATH
# wtdbg
export PATH=/home/u2052/jimmy/software/wtdbg-2.5_x64_linux:$PATH
# RepeatMasker
export PATH=/home/u2052/jimmy/software/RepeatMasker:$PATH
export PATH=/home/u2052/jimmy/software/RepeatMasker/util:$PATH
# RepeatModeler
export PATH=/home/u2052/jimmy/software/RepeatModeler-2.0.1:$PATH
# maker
export PATH=/home/u2052/jimmy/software/maker/bin:$PATH

```

#### 01.Braker训练物种的参数模型`Braker`
- 输入：
    - RNA-Seq.fastq（RNA-Seq数据）（建议是二代）
    - 组装的基因组序列（genome.fasta）
- 过程：
    - 参考基因组构建`hisat2-build`
    - 比对到比对到参考基因组`hisat2`
    - 训练模型`braker.pl`
- 输出：参考基因模型，放在augustus/config下
```
$ cat run_braker.sh
#!/bin/bash
# gm_key 文件如果过期需要重新下载，这个地方踩了坑，注意修改后的名字后面没有64，之前报错就是这个原因
cp /home/u2052/jimmy/software/gm_et_linux_64/gmes_petap/gm_key_64 ~/.gm_key

# 参数
species=sesame20200304      # 物种简写
genome=input/genome.fasta   # 基因组序列
threads=10                  # 线程数
rna_left=input/SRR7476717_1.fastq.gz
rna_right=input/SRR7476717_2.fastq.gz


# 为hisat2 构建参考序列
# hisat2-build $genome ./ 1>hisat2-build.log 2>&1

# 使用 hisat2 进行比对，若果有多个样本的数据，就运行多次比对，注意修改输入的 fastq 文件和输出的 bam 文件
# hisat2 --new-summary -p $threads -x ./ -1 $rna_left -2 $rna_right 2>hisat.log | samtools view -Sb - | samtools sort -o rnaseq.bam -

# 运行 braker
braker.pl --species=$species \
    --genome=$genome \
    --softmasking \
    --bam=rnaseq.bam \
    --cores=$threads \
    --AUGUSTUS_CONFIG_PATH=/home/u2052/jimmy/software/augustus-3.3.2/config \
    --BAMTOOLS_PATH=/home/u2052/jimmy/software/bamtools/bin \
    --SAMTOOLS_PATH=/home/u2052/jimmy/software/samtools \
    --GENEMARK_PATH=/home/u2052/jimmy/software/gm_et_linux_64/gmes_petap

```
#### 02.基因注释`Maker`
> 训练参数模型的时候用的转录组测序数据不需要拼接，直接将测序数据比对到组装好的基因组，用braker来训练，但是maker从头预测的时候需要拼接，三代直接全长转录本，一般用`Trinity`
>  Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G #左端的放一起，右端的放一起

- 输入：
    - 基因组序列（genome.fasta），组装好的
    - DNA层次的证据：训练的基因参数模型，在augustus config下
    - RNA层次序列（rna.fasta）
    - Protein层次证据（proteins.fasta），直接下载近缘物种的
- 输出：
    - gff/gtf
    - cds 序列
    - protein序列

- maker的使用有几个地方：第一个是安装，选择MPI，之后会差一个东西忘了，直接conda就行，然后重新 perl Build.PL 选择Y，之后可以直接用./build xxx这个命令安装maker需要的perl模块，但是大概率会报错，是因为有两个模块装不上一个是`DB_File``XML::LibXML::Reader`使用cpanm也装不上，然后从网上直接下载tar.gz，尝试安装，各种报错，很玄幻，其中DB_File的安装需要配置一下db.h，查了一下感觉是个数据库，之后就没有再折腾了，最后尝试重新perl Build.PL神奇的事情发生了，这一次没有报错，发现在maker目录下有一个bin大概率是安装成功了。
```
$ cat run_maker.sh
#!/bin/bash
# 将 maker 添加到 PATH
export PATH=/pub/software/maker-2.31.10/bin/:$PATH

# 生成配置文件
## maker -CTL

# 运行 maker
maker

# 合并成一个结果
fasta_merge -d genome.maker.output/genome_master_datastore_index.log  -o sesame
# -o：指定前缀，后缀自动生成
#sesame.all.maker.proteins.fasta ：蛋白序列
#sesame.all.maker.transcripts.fasta：cds序列
gff3_merge -d genome.maker.output/genome_master_datastore_index.log  -o sesame.all.maker.gff3
# -o：只有一个文件就直接指定
#sesame.all.maker.gff3  gff文件

# 修改基因、转录本ID
maker_map_ids --prefix SIN_ --justify 6 sesame.all.maker.gff3 > all.id.map
map_fasta_ids all.id.map sesame.all.maker.proteins.fasta
map_fasta_ids all.id.map sesame.all.maker.transcripts.fasta
map_gff_ids all.id.map sesame.all.maker.gff3

gffread -T -o sesame.all.maker.gtf sesame.all.maker.gff3 #转换为gtf
```
maker的参数极多，所以采用配置文件的方式，`maker -CTL`生成三个需要配置的文件
- 配置maker_opts.ctl
```
$ vim maker_opts.ctl
#-----Genome (these are always required)
#指定组装的genome位置
genome=../test_data/genome.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
#物种类型
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=../test_data/est.fasta #set of ESTs or assembled mRNA-seq in fasta format # 放自身
altest= #EST/cDNA sequence file in fasta format from an alternate organism # 近缘物种的
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=../test_data/protein.fasta  #protein sequence file in fasta format (i.e. from mutiple oransisms) # 多个物种的蛋白
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)#之前做过repeatmasker提供生成的物种重复序列的数据库，model_org删掉
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=/home/u2052/jimmy/workspace/02.Annotation/01.repeat/sesame-families.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/pub/software/maker-2.31.10/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction #下面的都是从头从基因预测软件，最好的就用augstus，更新快
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=sesame20200304 #Augustus gene prediction species model #输入构建的模型
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no #转录本对比到基因组预测
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no #蛋白比对到基因组预测
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no #小rna的预测
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs #核糖体rna要预测需要提供数据库
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no #这里有争议，频闭预测？还是非频闭预测？最好1

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI) #装过MPI选1

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage) #切分窗口的值越小越消耗资源（100K-1M）
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless) #最小的contig小于这个值的就不要算基因了

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=50 #require at least this many amino acids in predicted proteins #少于多少的氨基酸就不算基因，一般是50个
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no # 可变剪切参数，一般来说，对于非模式物种不要预测，除非测了大量的转录组数据
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no # 是否所有与的基因都是完整的，cds上翻译后没有启示密码子，终止密码子，错位等等，应该允许，因为这个是参考基因组
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```
- 配置maker_exe.ctl，maker_bopts.ctl
```
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/home/u2052/miniconda3/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/home/u2052/miniconda3/bin/blastn #location of NCBI+ blastn executable
blastx=/home/u2052/miniconda3/bin/blastx #location of NCBI+ blastx executable
tblastx=/home/u2052/miniconda3/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
RepeatMasker=/home/u2052/jimmy/software/RepeatMasker/RepeatMasker #location of RepeatMasker executable
exonerate=/home/u2052/miniconda3/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/usr/bin/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/home/u2052/jimmy/software/augustus-3.3.2/bin/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
tRNAscan-SE= #location of trnascan executable
snoscan= #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)
```
#### 03.BUSCO组装和注释的完整性评价
依赖数据库OrthrDB（直系同源数据库），选择物种的的时候尽量错开，每个层次都选一个物种，用busco的时候conda有点问题，一般来说直接用conda安装，因此用的是/pub/conda/busco
```
#!/bin/bash
# 使用不同数据库进行 BUSCO 分析
busco -m prot -i ./data/Sind.pep.fasta -o viridiplantae_odb10 -l database/viridiplantae_odb10 -c 4
        --------------------------------------------------
        |Results from dataset viridiplantae_odb10         |
        --------------------------------------------------
        |C:88.3%[S:85.3%,D:3.0%],F:3.0%,M:8.7%,n:430      |
        |380    Complete BUSCOs (C)                       |
        |367    Complete and single-copy BUSCOs (S)       |
        |13     Complete and duplicated BUSCOs (D)        |
        |13     Fragmented BUSCOs (F)                     |
        |37     Missing BUSCOs (M)                        |
        |430    Total BUSCO groups searched               |
        --------------------------------------------------

busco -m prot -i ./data/Sind.pep.fasta -o embryophyta_odb10 -l database/embryophyta_odb10 -c 4
busco -m prot -i ./data/Sind.pep.fasta -o eudicotyledons_odb10 -l database/eudicotyledons_odb10 -c 4

# 统计、比较、画图
mkdir BUSCO_summaries
cp */short_summary*.txt BUSCO_summaries/

generate_plot.py -wd BUSCO_summaries/ #这里可以去看一看是如何用python调用的R
```


### 基因功能注释： 
 > 通过基因注释的到了几万条蛋白的序列，就需要知道这些蛋白的功能是什么，如果只是少量的蛋白序列直接可以使用NCBI，Uniprot等在线数据库即可完成功能注释

- 同源注释：
    - 基于序列相似性（动态规划，滑动窗口，比对打分，前提假设，所 以位置的重要程度是一样的）：
        - 软件：`blast`,`Diamond`(大概14年，速度号称是blast的2W倍)
        - 数据库：Uniprot，Nr(再提一嘴，NCBI上非冗余蛋白数据库)
    - 基于结构域（计算了保守结构域的重要性）：
        - 软件：HMMER，Interproscan
        - 数据库：Pfam（针对结构域进行基因家族划分的基因家族数据库），interpro（他们存储了大量的重要结构域）

==> 通过同源比对（我们的哪一些序列跟数据库的哪些序列最为相似，为了知道这些序列的具体功能），在数据库中会对序列进行功能分类
- 功能分类：
    - GO：MF（分子功能）CC（细胞组分）BP（生物学过程） 
    - KEGG Pathway：研究很清楚的生物学过程（目前大概1000多个）
    - COG/KOG：基因家族分类，基于所有物种都是同源的，那么所有基因最早也有共同祖先分化而来，就形成了基因家族，这属于NCBI的体系，这是一个非常粗糙的基因家族，用A-Z表示所有的基因家族
基因功能注释方法EggNOG和之前转录组的一样：http://eggnog5.embl.de/#/app/home
```
# 对注释的文件进行整理，绘图，可以直接用一下张老师的emapperx.R这个程序，有时间可以研究一下写R包
$ cat emapperx.R
#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)

# read parameters
p <- arg_parser("make OrgDB from emapper")
p <- add_argument(p, "emapper_anno", help="emapper annotation result", type="character")
p <- add_argument(p, "proteins", help="proteins in fasta format", type="character")

argv <- parse_args(p)


library(tidyverse, quietly = TRUE)
library(formattable, quietly = TRUE)
library(AnnotationForge, quietly = TRUE)
library(seqinr, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)


# read emapper result
read_emapper <- function(emapper_anno) {
  emapper <- read_delim(emapper_anno, 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        comment = "#", trim_ws = TRUE) %>%
    dplyr::select(GID = X1, 
                  Gene_Symbol = X6, 
                  GO = X7, 
                  KO = X9, 
                  Pathway = X10, 
                  COG = X21, 
                  Gene_Name = X22)
  return(emapper)
}

anno_stat <- function(){
  
  # GO statistics and plot
  go_anno = length(unique(gene2go$GID))
  go_bp <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "BP",
                   level    = 2,
                   readable = FALSE)
  
  go_bp <- as.data.frame(go_bp)
  go_bp$GO_Class <- "Biological Process"
  
  go_cc <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "CC",
                   level    = 2,
                   readable = FALSE)
  
  go_cc <- as.data.frame(go_cc)
  go_cc$GO_Class <- "Cellular Component"
  
  go_mf <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "MF",
                   level    = 2,
                   readable = FALSE)
  go_mf <- as.data.frame(go_mf)
  go_mf$GO_Class <- "Molecular Function"
  
  go_all <- rbind(go_bp, go_cc, go_mf)
  write.table(go_all, "go.txt", sep = "\t", quote = F)
  p <- ggplot(go_all) + 
    geom_bar(aes(x = Description, 
                 y = Count,
                 fill = GO_Class),
             stat = "identity") + facet_wrap(~GO_Class, scales = "free_x") + 
    labs(title = "GO function classification", y = "Number of genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  ggsave("go.pdf", p, width = 20, height = 7)
  
  
  # Pathway statistics and plot ---------------------------------------------
  gene2pathway <- dplyr::select(emapper, GID, Pathway) %>%
    separate_rows(Pathway, sep = ',', convert = F) %>%
    filter(str_detect(Pathway, 'ko'))
  
  load(file = paste(script_dir, "kegg_info.RData", sep = "/"))
  
  gene2pathway <- gene2pathway %>%
    left_join(pathway2name, by = "Pathway") %>%
    dplyr::select(GID, Pathway, Pathway_Name, Pathway_Class, Pathway_Subclass) %>%
    distinct() %>%
    na.omit()
  
  pathway_anno = length(unique(gene2pathway$GID))
  
  pathway_stat <- dplyr::select(gene2pathway, GID, Pathway_Class, Pathway_Subclass) %>% 
    distinct() %>% 
    group_by(Pathway_Class, Pathway_Subclass) %>%
    summarise(Count = n(), Percentage = percent(n()/pathway_anno))
  
  pathway_stat$Pathway_Subclass <- ordered(pathway_stat$Pathway_Subclass, levels = pathway_stat$Pathway_Subclass) 
  
  p <- ggplot(pathway_stat, aes(x = Pathway_Subclass, y = Percentage)) +
    geom_bar(aes(fill = Pathway_Class), stat = 'identity') +
    geom_text(aes(label = Count), nudge_y = 0.005) +
    scale_y_continuous(labels=percent) + 
    labs(y = "Percent of genes(%)", x ="", fill = "Class") +
    coord_flip() +
    theme_classic()
  
  ggsave("pathway.pdf", p, width = 20, height = 7)
  write.table(gene2pathway, file = "pathway.txt", sep = "\t", quote = F)
  write.table(pathway2name, file = 'pathway_name.txt', sep = '\t', quote = F)
  write.table(pathway_stat, file = "pathway_stat.txt", sep = "\t", quote = F, row.names = F)
  
  
  # COG statistics and plot -------------------------------------------------
  library(readr)
  cog_funclass <- read_delim(paste(script_dir, "cog_funclass.tab", sep = "/"), 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
  
  insert_comma <- function(x){
    str_c(x, sep = '', collapse = ',')
  }
  
  gene2cog <- dplyr::select(emapper, GID, COG) %>%
    filter(!is.na(COG)) %>%
    mutate(COG = sapply(str_split(COG, ''), insert_comma)) %>%
    separate_rows(COG, sep = ',', convert = F) %>%
    left_join(cog_funclass, by = c('COG' = 'COG'))
  
  cog_anno = length(unique(gene2cog$GID))
  
  write.table(gene2cog, file = "cog.txt", sep = "\t", quote = F, row.names = F)
  
  p <- ggplot(data = gene2cog) + 
    geom_bar(aes(x = COG, 
                 fill = COG_Name)) +
    labs(title = "COG/KOG Function Classification ", 
         x = "",
         y = "Number of genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size=unit(1,"line"),
          legend.text = element_text(size = 7.5)) +
    guides(fill=guide_legend(ncol=1))
  ggsave("cog.pdf", p, width = 16, height = 7)
  
  # number and percentage ---------------------------------------------------
  anno_stat <- tibble(
    database = c("EggNOG", "GO", "COG/KOG", "KEGG Pathway"),
    number = comma(c(eggnog_anno, go_anno, cog_anno, pathway_anno), digits = 0),
    percentage = percent(c(eggnog_anno, go_anno, cog_anno, pathway_anno)/total_gene)
  )
  
  write.table(anno_stat, "anno_stat.txt", quote = F, row.names = F, sep = "\t")
}


# set script dir
script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

# get total gene number
all_gene <- getName.list(read.fasta(file = argv$proteins, 
                                    seqtype = 'AA'))
total_gene = length(all_gene)

emapper <- read_emapper(argv$emapper_anno)

# gene name
gene_info <- dplyr::select(emapper,  GID, Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name))
eggnog_anno = length(gene_info$GID)


# gene to gene ontology
gene2go <- dplyr::select(emapper, GID, GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  filter(!is.na(GO)) %>%
  mutate(EVIDENCE = 'IEA') 
  
# make org package
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               maintainer='zhangsan <zhangsan@genek.tv>',
               author='zhangsan',
               outputDir="./",
               tax_id=0000,
               genus='M',
               species='y',
               goTable="go",
               version="1.0")
  
pkgbuild::build('.//org.My.eg.db', dest_path = ".")

dir.create('R_Library', recursive = T)
install.packages('org.My.eg.db_1.0.tar.gz',
                   repos = NULL, #从本地安装
                   lib = 'R_Library') # 安装文件夹
library(org.My.eg.db, lib = 'R_Library')
anno_stat()

$ cat install.R
if (!requireNamespace("jsonlite", quietly = TRUE)) {install.packages("jsonlite")}
if (!requireNamespace("tidyverse", quietly = TRUE)) {install.packages("tidyverse")}
if (!requireNamespace("argparser", quietly = TRUE)) {install.packages("argparser")}
if (!requireNamespace("seqinr", quietly = TRUE)) {install.packages("seqinr")}
if (!requireNamespace("formattable", quietly = TRUE)) {install.packages("formattable")}

if (!requireNamespace("AnnotationForge", quietly = TRUE)) {BiocManager::install("AnnotationForge")}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {BiocManager::install("clusterProfiler")}
if (!requireNamespace("pathview", quietly = TRUE)) {BiocManager::install("pathview")}

$ Rscript emcp/emapperx.R query_seqs.fa.emapper.annotations Sind.pep.fasta
```
## => STEP3_比较基因组学：
- 什么是比较基因组学分析？（通过比较近缘物种的基因组（物种之间的差异，重测序是个体之间的差异））
    - 物种层次：这些物种的进化关系（单拷贝基因，从基因家族来）如何，是什么时候分化的（进化树构建与分歧时间估计）
    - 基因组层次：染色体片段是如何重组，复制，进化？（共线性与WGD）
    - 基因家族层次：基因家族有哪些差异（有和无，多和少）？（基因家族鉴定，特异，收缩，扩张）
    - 基因层次：基因结构，进化速率（正选择）
### 数据准备：
- 物种选择：
    - A.近缘物种比较：马铃薯，番茄，狸藻
    - B.生物特征比较
    - C.外群：葡萄
- 数据格式：
    - 蛋白序列：基因家族聚类，进化树
    - 基因注释：共线性
    - CDS序列：用于正选择
- 数据下载：
    - Phytozome：https://phytozome-next.jgi.doe.gov/
    - Ensembl：ftp://ftp.ensembl.org/pub/ ftp://ftp.ensemblgenomes.org/pub
    - NCBI
    - 物种数据库
- 数据要求：
    - 只保留最长转录本
    - CDS和PEP统一命名
    - fasta文件每行固定长度

**不演示代码了，总之要保证所有文件的id是一一对应的，可能会遇到下面一种情况**
```
>AAA    BBB    CCC
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>DDD    EEE    FFF
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# 需求：将含有>的行中CCC这一列替换到AAA这一列的位置，其他行保持不变
# 张老师推荐用
$ map_fasta_ids id.map.txt file.fasta #maker中的一个函数 #maker中的一个函数
$ cat id.map.txt # id替换表
AAA    CCC
DDD    FFF

# 觉得可以用awk写函数完成，试了一下，不知道awk怎么锁定以>开头的行，NR显示的是行号，从1开始计数
# 因此需要将原pep.fasta转化为一行用到
$ fasta_formatter -i file.pep.fasta -o sub.pep.fasta -w 1000000 #转化为一行
$ awk '{tmp=NR%2;if(tmp == 1)print ">"$3;else print $0}' sub.pep.fasta | less #也可以实现效果，python当然也可以不过效率肯定没有awk高
```
### 基因家族聚类：
- OrthoVenn2基于著名的OrthoMCL，直接网页版就很好用
- OrthoFinder近几年发布，效果应该比前者好L

>  物种的所有pep.fasta文件，diamond会任意的两两进行比对（msa），得到哪些序列是相似的，接着通过MCL的方法进行聚类，得到Ortholog group（基因家族），MCL同时构建genetree进化树（基因树，对每一个基因家族的所有基因进行构建进化树），最后通过构建的所有的基因家族的基因树再构建物种树（sptree）

> genetree的构建：一种是基于MSA（多序列比对，把所有序列一一对齐（左右端对齐，对不齐的填空）进行全局比对）的；一种是不基于MSA的，速度快，但不推荐

> spstree的构建：一种是基于genetree的；另一种是基于连在一起的序列（单拷贝基因连在一起）

> 树文件的格式：`file.nwk`可以通过MAGA打开

> 基因家族的注释：会冲突，可能基因注释的功能不一样，很难定义一个基因家族的功能

> 关于基因树和物种树的若干疑问：大部分基因家族的关系和物种树是一样的，因为基因的分开是伴随着物种的分开而分开，有个别基因树与物种树不一致，有可能是错误（序列太短），有可能是正选择，某个物种的这个基因受到了正选择，演化速度变快，与其他物种就渐行渐远了，还有一种可能，是基因的水平转移，比如植物和真菌长期共生的过程过，基因发生了水平转移（一般情况下基因是垂直转移，物种的分化繁殖传递到下一代），这种现象研究的不透彻
```
# 不懂为什么\后面不能这样注释，真正跑的时候写成一行
$ orthofinder -f data \ #输入，所有物种的pep.fasta文件都放在data下面，注意命名，file.fasta别太长
    -S diamond \  # 比对软件用diamond
    -M msa \ # 多序列比对
    -T fasttree \
    -t 20
```
### 物种层次：
- 进化关系：input：单拷贝基因；output：物种树（sp.nwk）
- 分化时间：
> 构建一个物种树：基于所有物种共有的单拷贝基因家族，两种方法：1.连接（越长的基因对建树的影响越大）：将多个物种公有的按拷贝基因家族各自分别连接起来，进行多序列比对，然后对每个物种进行建树；2.合并（每条基因平等对待），对每个物种的每个公有的单拷贝基因家族进行分别构树，接着汇总整合得到一致树。（不论怎么建树，都是先多序列比对然后用ML或者BI进行建树）
#### 01.进化树构建：
（基于基因家族聚类的结果来做，涉及到太多算法，感觉了解一下即可）
**进化树**：系统发育或系统发育树（phylogenetic tree），也叫进化树，是物种间，基因间，群体间，乃至个体间谱系关系的一种表现形式。

- **有根树和无根树**：根据是否指定了根节点，系统发育树可分为有根树和无根树
    - 无根树：没有指定祖先节点：只能看出各个节点的拓扑结构和相对距离
    - 有根树：反映了树上物种或基因的时间顺序，一般采用**外群定根法**，建树时引入亲缘关系较远的物种作为外群来设定根

- **树的格式**：1.Newick format；2.The New Hampshire X Format(NHX)；3.Nexus format
- **树的构建方法**：
    - 基于距离/分歧度：（两两比对算差异）
        - 非加权算术平均对群法UPGMA
        - 邻接法Neighbor-joining 
    - 基于特征/性状：具体的碱基序列和氨基酸序列
        - 最大简约法：最小变化树（祖先状态最小化）
        - 最大似然法：所有枝长和模型参数最优化
        - 贝叶斯推断：基于后验概率（用枝长和后验概率联合计算）
- **自展检验**：用来检验所计算的进化树分支可信度

**三种建树的方法**：iTOL，Orthofinder(基于连接好的supergene)，Astral（基于合并）
```
# raxml构树：最大似然法
# 安装直接conda
$ muscle   -in all_FAD.fasta > all_FAD.muscle.fasta #首先进行多序列比对
# 如果数据量太大，多序列比对可以用 mafft
$ raxmlHPC-PTHREADS -T 15 -m PROTGAMMAJTT -f a -p 123 -x 123 -# 100 -n out -s all_FAD.muscle.fasta  1>tree.log 2>tree.err # 使用raxmlHPC这个是多线程的
# -m 指定模型，这个看情况要用的时候研究一下，iq-tree可以找一找模型
# -f a：快速的构树
# -p：随机数种子，抽取
# -x：随机数种子，检验
# -#：做100次检验
# -n：指定输出文件的后缀名
# -s：指定输入的多比对信息

# 物种树的构建：
#Orthofinder结果目录中多序列比对目录下的每个物种的单拷贝基因家族连接在一起的
$grep -A 2 '>' SpeciesTreeAlignment.fa 
>Sind
MWRNSLSLSSKLASPNSNPSLLLLLHHGLLARSNSTS----PESTSAADPV-------SWAQSGTPRPSASDFPSLAAAA
AVRAKEVADLAKHYGRCYWELSKARLSMLVVATSGTGYILGSGSSIDYLGLCCTCAGTMMVAASANSLNQVFEVKNDAKM
--
>Slyc
MWRNSVGFSSKLLTSKSSCSRLPLLNHGVASRSTSTAAGSGPESTEAAGTVKVGFSSTEWVRTGSPTSSAVDLSSFAAAS
VLKGRDIVDMAKHYGRCYAELSKARLSMLVVATSGTGYILGSGSAIDYMGLCCTCAGTMMVAASASTLNQVFEVKNDAKM
--
>Stub
MWRNSVGFSSKLLTSKSSYSRLPLLNHGVASRSTSTSAGSGPESTEAAGTVKVGFSSTEWVRTGSPGSSAVDLSSFAAAS
VLKGRDIVDMAKHYGRCYAELSKARLSMLVVATSGTGYILGSGSAIDYMGLCCTCAGTMMVAASASTLNQVFEVKNDAKM
--
>Vvin
MWRNSQSLYSKLVFSRNPNPNPTFL---TSSSISSLDAIVRPFSIASDGTRTVG--------------STINTTSLSA--
----RDAVDLARHYGRCYCELSKARLSMLVVATSGTGFVLGSGNIIDFGGLFWTCAGTMMVAASANSLNQVFEINNDAKM
$ raxmlHPC-PTHREADS -T 14 -m PROTGAMMAJTT -f a -p 123 -x 123 -# 100 -n out -s SpeciesTreeAlignment.fa 1>tree.log 2>tree.err

# Astral，每个单拷贝基因家族单独建树文件
# 单拷贝基因家族的列表
$ cat ~/jimmy/workspace/04.GeneFamilyCluster/data/OrthoFinder/Results_Apr28/Orthogroups/Orthogroups_SingleCopyOrthologues.txt 
# 每个基因家族的树在Gene_Trees中，于是要做的就是，基于单拷贝基因家族列表，从Gene_Trees中提取出来
$ cat run_astral.sh
#!/bin/bash

genetree="/home/u2052/jimmy/workspace/04.GeneFamilyCluster/data/OrthoFinder/Results_Apr28/Gene_Trees"
sco="/home/u2052/jimmy/workspace/04.GeneFamilyCluster/data/OrthoFinder/Results_Apr28/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

cat $sco | while read lines
do
        cat $genetree/$lines\_tree.txt | awk '{print $0}' 
done > SingleCopyTrees

sed 's/_[a-zA-Z]\+.\{0,25\}:/:/g' SingleCopyTrees > Astral.input.trees #去掉基因名保留物种名和枝长
#sed 's/\([(,]\{1\}[A-Za-z]\+\)_[^:]\+/\1/g' SingleCopy.trees > Astral_input.trees #王莹老师的写法

```
#### 02.分歧时间估计：
- 基于分子钟假说：氨基酸或核苷酸替代速率在进化过程中近似地保持恒定：
    - 全局分子钟模型：序列间的期望距离随分化时间线性增加
    - 局部分子钟模型：对分支赋予不同的进化速率
### 基因组层次：
- 共线性：input：基因注释.gff，蛋白序列RAxML_bipartitionsBranchLabels.out
### 基因家族：
- 基因家族的鉴定：输入（蛋白质序列）；输出：基因家族 --单拷贝基因
- 有和无
- 多和少：扩张与收缩
### 基因层次：
- 基因结构
- 进化速率：正选择
