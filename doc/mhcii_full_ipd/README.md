# MHC-II IPD 序列清理说明

本目录保存从 IPD 原始蛋白 FASTA 重建的 MHC-II 链表，供
`convert_tcr_pmhc_full_sequences.py` 生成 TCR-pMHC-II 结构预测输入。

## `full_aa` 的定义

`full_aa` 不是包含信号肽和跨膜区的完整前体，而是用于可溶性复合物建模的
MHC-II 成熟结构域核心：

- alpha 链保留 alpha1/alpha2 结构域；
- beta 链保留 beta1/beta2 结构域；
- 去除 N 端信号肽；
- 去除 C 端连接肽、跨膜螺旋和胞质尾。

这一定义与本项目 MHC-I 表的处理方式一致，可避免 AlphaFold/复合物建模时把
信号肽或膜锚当成复合物中的游离结构。`full_aa_original` 和 `sequence` 保存 IPD
原始蛋白序列，便于追溯。

## 输出文件

| 文件 | 用途 |
| --- | --- |
| `mhcii_chain_qc_all.csv` | 所有支持物种、支持基因的原始记录及逐条 QC 结果 |
| `mhcii_chain_curated.csv` | 通过成熟链完整性检查的全部精确 IPD 记录 |
| `mhcii_model_input.csv` | 转换脚本默认读取的安全查询表，每个无歧义名称键一条序列 |
| `mhcii_two_field_conflicts.csv` | 同一二字段名称对应多条成熟蛋白序列的隔离记录 |
| `mhcii_cleaning_summary.json` | 阈值、参考链、原始文件 SHA256、计数和拒绝原因 |

当前原始文件共得到：

- 15,549 条支持范围内的 MHC-II 记录；
- 6,440 条完整且可用于建模的精确来源记录；
- 2,868 个无蛋白序列歧义的查询键；
- 3 个有歧义的非人源二字段键被排除出模型输入表。

计数以 `mhcii_cleaning_summary.json` 为准；更换 IPD 原始版本后会随重建结果变化。

## 重建方法

脚本直接读取：

```text
raw/IPD-IMGT-HLA_hla_prot.fasta
raw/IPD-MHC_MHC_prot.fasta
```

不使用旧版 `mhcii_model_input.csv` 作为边界来源。六类成熟结构域参考如下，坐标均
为原始前体序列的 0-based、右端不包含坐标：

| 家族 | 链 | 参考等位基因 | `[start, end)` | 长度 |
| --- | --- | --- | ---: | ---: |
| DR | alpha | `HLA-DRA*01:01:01:01` | `[25, 207)` | 182 |
| DR | beta | `HLA-DRB1*01:01:01:01` | `[29, 220)` | 191 |
| DQ | alpha | `HLA-DQA1*01:01:01:01` | `[23, 207)` | 184 |
| DQ | beta | `HLA-DQB1*05:01:01:01` | `[32, 223)` | 191 |
| DP | alpha | `HLA-DPA1*01:03:01:01` | `[31, 211)` | 180 |
| DP | beta | `HLA-DPB1*01:01:01:01` | `[29, 217)` | 188 |

每条候选序列与相同家族、相同链型的参考成熟链做局部蛋白比对，再转移成熟链
边界。端点只允许极小的参考端未对齐，并且只有候选原始序列确实存在对应残基时
才补回；不会对缺失片段进行填充或推断。

记录必须同时满足：

- 仅包含标准 20 种氨基酸；
- 不带 `N/L/S/C/A/Q` 表达状态后缀；
- 参考覆盖率至少 0.97；
- 对齐区一致性至少 0.72；
- 参考 N 端未对齐不超过 2 aa，C 端不超过 3 aa；
- 成熟链长度与对应参考相差不超过 5 aa；
- 任意 19-aa 窗口最多含 13 个疏水残基，防止跨膜段残留。

`qc_reason` 会列出每条拒绝记录触发的全部原因。

## 名称消歧

转换脚本允许完整等位基因名称回退到二字段名称。为避免旧表中“按文件顺序取
第一条”的隐式行为，重建时会按二字段键检查 `full_aa`：

- 同一键的所有精确记录具有相同蛋白序列时，在 `mhcii_model_input.csv` 中合并为
  一条；`source_name` 指向实际选用的 IPD 记录；
- 同一键对应多条蛋白序列时，整组不进入模型输入表，并写入
  `mhcii_two_field_conflicts.csv`；
- 精确来源记录仍保留在 `mhcii_chain_curated.csv`，不会丢失审计信息。

当前隔离的键为：

```text
Mamu-DRB1*07:03
Patr-DPB1*01:13
Patr-DRB1*03:07
```

## 运行

在项目根目录执行：

```bash
python -m pip install pandas biopython
python doc/mhcii_full_ipd/rebuild_mhcii_tables.py
```

脚本会重建上述 CSV/JSON，并在写出前执行内部一致性检查。也可以先写到临时目录：

```bash
python doc/mhcii_full_ipd/rebuild_mhcii_tables.py \
  --output-dir /tmp/mhcii_check
```

## 依据与限制

原始序列来自 [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) 和
[IPD-MHC](https://www.ebi.ac.uk/ipd/mhc/download/)。成熟前体、胞外区、跨膜区和
连接肽的注释可在 UniProt 的
[HLA-DRA](https://www.uniprot.org/uniprotkb/P01903/entry)、
[HLA-DRB1](https://www.uniprot.org/uniprotkb/P01911/entry)、
[HLA-DQA1](https://www.uniprot.org/uniprotkb/P01909/entry)、
[HLA-DQB1](https://www.uniprot.org/uniprotkb/P01920/entry)、
[HLA-DPA1](https://www.uniprot.org/uniprotkb/P20036/entry) 和
[HLA-DPB1](https://www.uniprot.org/uniprotkb/P04440/entry) 条目中核对。

本表的 `full_aa` 特意终止于膜近端结构域末端，不包含 UniProt 注释的柔性连接肽，
因此应称为“结构预测用成熟核心链”，不能称为生物学意义上的完整前体或完整胞外
区。跨物种边界通过同家族人源参考转移；低于严格覆盖率/一致性阈值的远缘或部分
序列会被保守地拒绝。
