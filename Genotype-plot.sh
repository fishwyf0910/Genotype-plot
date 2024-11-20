conda activate vcftools
# 提染色体
vcftools --gzvcf input.vcf --chr n --recode --recode-INFO-all --stdout | bgzip > output.vcf.gz

# 提样本
bcftools view -s 02R0004,03R0016,02R0010,05R0020,03R0028,02R0018,03R0008,05R0006,03R0043,05R0032,03R0027,03R0001,XG-1,XG-2,XG-3,HB-1,HB-2,HB-3,HB-4,HB-5,HB-6,OJ-1,OJ-2,OJ-3 /data01/wangyf/project2/CyprinusCarpio/15.pop/0.vcfdata/final-raw.indel5.biSNP.QUAL30.QD3.FS20.MQ55.SOR3.MQRS-5.RPRS-5.PASS.GQ10.popmiss.maxmiss0.15.AF0.05.10-3ClusterFilter.vcf.gz -Ov  > haplo/chrA08/abtb1/2pop.vcf

bgzip haplo/chrA08/abtb1/2pop.vcf && tabix -p vcf haplo/chrA08/abtb1/2pop.vcf.gz

vcftools --gzvcf 2pop.vcf.gz --positions specific_position.txt --recode --out specific_position.vcf

# 提基因且受选择的区域
g=abtb1; mkdir ../$g; cd ../$g
# abtb1
bcftools filter /data01/wangyf/project2/CyprinusCarpio/15.pop/18.dtest/irtyshsouth/haplo/chrA08/2pop.vcf.gz --regions NC_056579.1:4157263-4168488 > $g.vcf
# mrpl40:       NC_056579.1:2443916-2443936
# ubox5:        NC_056579.1:988479-988520
# naif1:        NC_056579.1:2476800-2476800
# asic1a:       NC_056579.1:4125581-4127775
# grip2a:       NC_056579.1:4194175-4222364
# sema6a:       NC_056579.1:601712-603311
# plat: NC_056579.1:2347186-2351504
# ikbkb:        NC_056579.1:2378830-2397730
# gle1: NC_056579.1:2454112-2458171
# lrrc2:        NC_056579.1:869228-871827
# prr16:        NC_056579.1:666221-666382
# snx2: NC_056579.1:697739-716100

# vcf格式转plink格式
plink --vcf $g.vcf --recode --out $g --const-fid --allow-extra-chr
plink --allow-extra-chr --file $g --noweb --make-bed --out $g

# 用plink将基因型转成0、1、2格式
plink --bfile $g --allow-extra-chr --recodeA --out $g

# 删除第1行、第1、3、4、5、6列
sed '1d' $g.raw | awk '{ $1=""; $3=""; $4=""; $5=""; $6=""; print $0 }' | sed 's/^ *//; s/  */ /g' > $g.plot

# 添加表头
{ echo "header $(head -1 $g.plot | awk '{for(i=2;i<=NF;i++) printf("t%d%s", i-1, (i<NF?" ":""))}')"; cat  $g.plot ; } > $g.r.plot
sed -i "s/ /\t/g" $g.r.plot
# 手动删掉header

conda activate r4.3
R
library(pheatmap)

# 读取数据
test <- read.table("138snp.allsample.r.plot", header = TRUE, row.names = 1)

# 读取name.list
name_list <- readLines("name.list")
test_ordered <- test[name_list, , drop = FALSE]

# 绘制热图
pdf("ik.5site.specific_position.pdf", width = 3, height = 6)  # 设置宽度为8， 高度为6

pheatmap(test_ordered, 
         color = colorRampPalette(c("white", "#8bb5d1", "#cb5c5b"))(3), 
         cluster_col = FALSE, 
         cluster_row = FALSE, 
         show_colnames = FALSE)

dev.off()

# 绘制 EHH 图
# 生成两个群体的vcf文件
cd /data01/wangyf/project2/CyprinusCarpio/15.pop/18.dtest/irtyshsouth/haplo/chrA08/abtb1

bcftools view -s XG-1,XG-2,XG-3,HB-1,HB-2,HB-3,HB-4,HB-5,HB-6,OJ-1,OJ-2,OJ-3 /data01/wangyf/project2/CyprinusCarpio/15.pop/0.vcfdata/final-raw.indel5.biSNP.QUAL30.QD3.FS20.MQ55.SOR3.MQRS-5.RPRS-5.PASS.GQ10.popmiss.maxmiss0.15.AF0.05.10-3ClusterFilter.vcf.gz -Ov  > south.vcf

bgzip south.vcf && tabix -p vcf south.vcf.gz

vcftools --gzvcf south.vcf.gz --chr NC_056579.1 --recode --recode-INFO-all --stdout | bgzip > south.chrA08.vcf.gz

# 先phasing
conda activate beagle
java -jar /data01/wangyf/software/beagle.06Aug24.a91.jar gt=south.chrA08.vcf.gz out=south.chrA08.phasing

# 把snp名称补上
zcat south.chrA08.phasing.vcf.gz | awk 'BEGIN { OFS="\t" } /^#/ { print; next } { $3 = $2; print }' > south.chrA08.phasing.snp.vcf

# 去RStudio
setwd("E:/Rworkspace/EHH/")
library(rehh)
library(R.utils)
library(ggplot2)

hh <- data2haplohh(hap_file = "south.chrA08.phasing.snp.vcf",polarize_vcf = FALSE,vcf_reader = "data.table")

# abtb1
res <- calc_ehh(hh,mrk = "4167004",limhaplo = 2,include_zero_values=TRUE,include_nhaplo = TRUE)
  plot(res,main = NULL,xlab = "A08 (MB)",ylab="EHH",xlim=c(4162000,4172000),col =c("#ffb3c1","#e53935"),lwd=1,cex.axis=0.8,cex.lab=1)

pdf("rehh.south.A08.abtb1.pdf",width=4,height=4)
  plot(res,main = NULL,xlab = "A08 (MB)",ylab="EHH",xlim=c(4162000,4172000),col =c("#ffb3c1","#e53935"),lwd=1,cex.axis=0.8,cex.lab=1)
# plot(res,main = NULL,xlab = "A08 (MB)",ylab="EHH",xlim=c(4162000,4172000),col =c("#ade8f4","#2AA5DC"),lwd=1,cex.axis=0.8,cex.lab=1)
dev.off()

#######################################################################################################################################
# 得到hapbin的input文件
conda activate vcftools
vcftools --vcf south.1mb.vcf.phasing.snp.vcf --IMPUTE --out south.1mb.vcf.phasing.snp

plink --vcf south.1mb.vcf.phasing.snp.vcf --recode --out south.1mb.vcf.phasing.snp.impute --const-fid --allow-extra-chr

/data01/wangyf/software/hapbin/build/ehhbin --hap south.1mb.vcf.phasing.snp.impute.hap --map south.1mb.vcf.phasing.snp.impute.map --locus 4167012 > south.1mb.vcf.phasing.snp.impute.ehh

# ehhbin输出五列。前三个是位点的ID及其遗传和物理位置。接下来是两列，对应于该位点每个等位基因的EHH
{ echo -e "pos\tref\talt";   awk '!/IHH|iHS|MAF/ {print $1 "\t" $4 "\t" $5}' south.1mb.vcf.phasing.snp.impute.ehh; } > south.1mb.vcf.phasing.snp.impute.ehh.plot

conda activate r4.3
R
library(ggplot2)

data <-read.table('south.1mb.vcf.phasing.snp.impute.ehh.plot',header=TRUE)

p <- ggplot(data, aes(x = pos)) +
  geom_line(aes(y = ref, color = "Ref"), size = 1) +  # 第一条折线
  geom_line(aes(y = alt, color = "Alt"), size = 1) +  # 第二条折线
  labs(x = "Position",y = "EHH") +
  xlim(3600000, 4600000) + ylim(0, 1) +
  scale_color_manual(name = "Legend", values = c("Ref" = "#87CEFA", "Alt" = "#1E90FF")) +
# scale_color_manual(name = "Legend", values = c("Ref" = "#FFB6C1", "Alt" = "#DC143C")) +
  theme_minimal()

ggsave("south.1mb.vcf.phasing.snp.impute.ehh.plot.pdf", plot = p, device = "pdf", width = 7, height = 6)
