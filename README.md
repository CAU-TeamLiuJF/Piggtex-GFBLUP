# PigGTEx-GFBLUP
基于 PigGTEx QTL 注释的猪基因组预测（GFBLUP）分析。

> 注：以下以 lncQTL（data5 群体，性状 AGE/BF/TNB）为例；切换其他 QTL 类型或数据集时，
> 需同步修改各脚本顶部写死的 PigGTEx 目录、rds 编号与 traits。

## 脚本功能

| 脚本 | 功能 |
|---|---|
| `piggtex_d.R` | 将 PigGTEx 显著 QTL 与芯片 SNP 匹配（精准 + ±50 kb），按组织输出 QTL SNP 列表 |
| `qtl_piggtex_gs_50kb.R` | 用 QTL SNP + 背景 SNP 双随机效应 GFBLUP 做 5 折×2 重复交叉验证，计算预测准确性 cor 与 bias |
| `result_summary.R` | 汇总各 性状×组织 的交叉验证结果，输出 cor/bias 均值±SD |
| `submit_qtl.sh` | SLURM 数组作业批量提交（每任务 = 一个 性状×组织 组合） |

## 使用方法

```bash
# 1. 生成 QTL 注释（参数 1~5 对应不同芯片，见脚本注释）
Rscript piggtex_d.R 5

# 2. 提交交叉验证（先确认 submit_qtl.sh 中启用的 Rscript 行）
sbatch submit_qtl.sh

# 3. 汇总结果
Rscript result_summary.R
```

依次运行即可。
