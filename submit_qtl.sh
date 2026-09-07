#!/bin/bash
#SBATCH -J significant_vc_lncQTL_data5
#SBATCH -o log/%A_%a_significant_vc_lncQTL_data5.out  # 使用作业ID和任务ID命名输出文件
#SBATCH --array=1-102%5             # 创建102个任务，同时最多运行10个 68, 170   #1-11%11   -4%4
#SBATCH -N 1                          # 每个任务使用1个节点
#SBATCH --cpus-per-task=10             # 每个任务分配4个CPU核心
#SBATCH --mem=150G                      # 每个任务分配20GB内存
#SBATCH -p Cnode_all 
#SBATCH --exclude=cnode2041

set -e
module unload R
source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3
cpath=/public/home/liujf/workspace/shikp/TempWork/20250612_piggtex/code
# 获取当前任务ID
TASK_ID=$SLURM_ARRAY_TASK_ID

# 格式化任务ID为3位数字（保持与原始seq -w行为一致）
FORMATTED_ID=$(printf "%03d" $TASK_ID)

# 执行R脚本
#Rscript ${cpath}/qtl_piggtex_gs_50kb.R ${FORMATTED_ID}
#Rscript ${cpath}/piggtex_d.R ${FORMATTED_ID}
#Rscript ${cpath}/eqtl_piggtex_variance.R ${FORMATTED_ID}
#Rscript ${cpath}/eqtl_gblup_variance.R ${FORMATTED_ID}
#Rscript ${cpath}/gwas_gs.R ${FORMATTED_ID}
#Rscript ${cpath}/gwas_gs_p.R ${FORMATTED_ID}
Rscript ${cpath}/eqtl_variance_decomp.R ${FORMATTED_ID}
echo "Task ${FORMATTED_ID} completed!"
