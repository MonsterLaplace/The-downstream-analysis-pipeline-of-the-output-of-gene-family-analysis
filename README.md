## **The-downstream-analysis-pipeline-of-the-output-of-gene-family-analysis**

##**First install seqkit**

mamba install -c bioconda seqkit

##**Expansion gene family analysis** **扩张基因家族基因提取**

# 提取Gamma_change.tab第7列代表物种PRP的收缩的orthogroupsID
cat Gamma_change.tab | awk '{if($6>0) print $0}' | cut -f1,6 | awk '{print $1"\t"$2}' >> TAntelope.expanded
# 根据sample ID和编号提取sample分支的基因家族显著扩张或收缩的基因家族树（Gamma_asr.tre文件中默认以p<0.05为标准判断变化是否显著）
grep "TAntelope<5>\*" Gamma_asr.tre > TAntelope_significant_trees.tre
# 提取sample分支显著变化的OG IDs （默认以p<0.05为标准）
grep -E -o "OG[0-9]+" TAntelope_significant_trees.tre > TAntelope_significant.ogs
# 以p<0.01为标准提取所有显著扩张或收缩的orthogroupsID（根据情况调整，常用p<0.05或p<0.01）
awk '$2 <0.01 {print $1}' Gamma_family_results.txt >p0.01_significant.ogs
# 提取以p<0.01为标准判断显著性的sample分支基因家族显著变化的OG IDs
grep -f TAntelope_significant.ogs p0.01_significant.ogs > TAntelope_p0.01_significant.ogs
#提**取显著收缩的sample物种的orthogroupsID**
grep -f TAntelope_p0.01_significant.ogs TAntelope.expanded | cut -f1 > TAntelope.expanded.significant
#**提取显著收缩的基因列表，假设基因ID是PRP_的前缀**
grep -f TAntelope.expanded.significant Orthogroups.txt | sed "s/ /\n/g"|grep "TAntelope" |sort -k 1.3n | uniq > TAntelope.expanded.significant.genes
#**提取显著收缩的基因序列**
seqkit grep -f TAntelope.expanded.significant.genes /mnt/z/xb/orthofinder/TAntelope.faa >TAntelope.expanded.significant.pep.fa

##**Contraction gene family analysis****收缩基因家族基因提取**

# 提取Gamma_change.tab第7列代表物种PRP的收缩的orthogroupsID
cat Gamma_change.tab |cut -f1,8|grep "-" >TAntelope.contracted
# 根据sample ID和编号提取sample分支的基因家族显著扩张或收缩的基因家族树（Gamma_asr.tre文件中默认以p<0.05为标准判断变化是否显著）
grep "TAntelope<5>\*" Gamma_asr.tre > TAntelope_significant_trees.tre
# 提取sample分支显著变化的OG IDs （默认以p<0.05为标准）
grep -E -o "OG[0-9]+" TAntelope_significant_trees.tre > TAntelope_significant.ogs
# 以p<0.01为标准提取所有显著扩张或收缩的orthogroupsID（根据情况调整，常用p<0.05或p<0.01）
awk '$2 <0.01 {print $1}' Gamma_family_results.txt >p0.01_significant.ogs
# 提取以p<0.01为标准判断显著性的sample分支基因家族显著变化的OG IDs
grep -f TAntelope_significant.ogs p0.01_significant.ogs > TAntelope_p0.01_significant.ogs
#**提取显著收缩的sample物种的orthogroupsID**
grep -f TAntelope_p0.01_significant.ogs TAntelope.contracted | cut -f1 > TAntelope.contracted.significant
#**提取显著收缩的基因列表，假设基因ID是PRP_的前缀**
grep -f TAntelope.contracted.significant Orthogroups.txt | sed "s/ /\n/g"|grep "TAntelope" |sort -k 1.3n | uniq > TAntelope.contracted.significant.genes
#**提取显著收缩的基因序列**
seqkit grep -f TAntelope.contracted.significant.genes /mnt/z/xb/orthofinder/TAntelope.faa >TAntelope.contracted.significant.pep.fa
