# TCR-pMHC 全链氨基酸序列转换

本说明对应脚本：

```text
convert_tcr_pmhc_full_sequences.py
```

脚本读取包含 TCR 基因、CDR3、MHC 和 peptide 信息的 CSV，将每一行转换为五条用于复合物结构建模的氨基酸链：

| 输出列 | 内容 |
| --- | --- |
| `A` | TCR alpha 链 |
| `B` | TCR beta 链 |
| `C` | MHC-I 重链，或 MHC-II alpha 链 |
| `D` | MHC-I 的 beta-2-microglobulin（B2M），或 MHC-II beta 链 |
| `E` | peptide |

输出 CSV 只包含上述五列，并严格保持 `A,B,C,D,E` 的顺序。

## 1. 运行环境

建议使用 Python 3.10 或更高版本，并安装 `pandas`：

```bash
python -m pip install pandas
```

项目中已经包含默认使用的参考序列表：

```text
doc/tcr_gene_aa_human_TRA_TRB_gene.csv
doc/mhci_full_ipd/mhci_full_ipd_curated_model_input_species_matched_b2m.csv
doc/mhcii_full_ipd/mhcii_model_input.csv
```

默认 TCR 参考表是人源 `TRA/TRB` 基因表。因此，除非通过参数提供其他参考表，否则 TCR 基因应为人源 `TRAV/TRAJ/TRBV/TRBJ`。

## 2. 快速开始

在项目根目录运行：

```bash
python convert_tcr_pmhc_full_sequences.py input.csv -o output.csv
```

如果不指定 `-o`，脚本会在输入文件同目录生成：

```text
<输入文件名>_full_sequences.csv
```

例如：

```bash
python convert_tcr_pmhc_full_sequences.py data/my_tcr.csv
```

默认输出为：

```text
data/my_tcr_full_sequences.csv
```

查看全部命令行参数：

```bash
python convert_tcr_pmhc_full_sequences.py --help
```

## 3. 输入 CSV

每一行必须提供以下八类信息。列名不区分大小写，并兼容项目中常见的别名；同一类信息只能保留一个列名。

| 信息 | 推荐列名 | 兼容列名 |
| --- | --- | --- |
| MHC | `MHC` | `mhc`, `HLA`, `hla`, `mhc_name` |
| peptide | `peptide` | `Peptide`, `PEPTIDE`, `epitope`, `Epitope` |
| TCR alpha CDR3 | `CDR3A` | `cdr3a`, `CDR3a`, `cdr3a_aa` |
| TCR beta CDR3 | `CDR3B` | `cdr3b`, `CDR3b`, `cdr3b_aa` |
| alpha V gene | `TRAV` | `trav`, `av_gene` |
| alpha J gene | `TRAJ` | `traj`, `aj_gene` |
| beta V gene | `TRBV` | `trbv`, `bv_gene`, `v_gene` |
| beta J gene | `TRBJ` | `trbj`, `bj_gene`, `j_gene` |

推荐的输入表头：

```csv
MHC,peptide,CDR3A,CDR3B,TRAV,TRAJ,TRBV,TRBJ
```

MHC-I 输入示例：

```csv
MHC,peptide,CDR3A,CDR3B,TRAV,TRAJ,TRBV,TRBJ
HLA-B07:02,HPNGYKSLSTL,CAYRSAGSGGSNYKLTF,CASSLMTNQPQHF,TRAV38-2/DV8,TRAJ53,TRBV27,TRBJ1-5
HLA-B*07:02,HPNGYKSLSTL,CAYRSATGYALNF,CASSLYGSVQNEQFF,TRAV38-2/DV8,TRAJ41,TRBV27,TRBJ2-1
```

输入要求：

- `CDR3A` 和 `CDR3B` 必须是氨基酸序列，并从保守的 `C` 开始。
- TCR 基因必须分别属于 `TRAV`、`TRAJ`、`TRBV` 和 `TRBJ` 家族。
- 氨基酸序列会转换为大写，并移除空格、`.` 和 `-`。
- 允许的残基字符为标准 20 种氨基酸和 `X`；其他字符会被拒绝。
- 基因名可带等位基因后缀，例如 `TRBV7-9*01`。精确名称不存在时，会尝试去掉 `*01` 后查询。
- 输入文件可以包含额外列，但额外列不会写入五列输出文件。

## 4. MHC-I 输入规则

MHC-I 的 `MHC` 单元格只需提供一个等位基因，例如：

```text
HLA-A02:01
HLA-A*02:01
HLA-B07:02
HLA-C*03:04
```

脚本会自动输出：

- `C`：对应 MHC-I 重链；
- `D`：参考表中物种匹配的 B2M 序列。

MHC 名称优先按完整等位基因精确匹配；没有精确记录时，会回退到 two-field 名称，例如将 `HLA-A*02:01:01:01` 和 `HLA-A*02:01` 归入 `HLA-A*02:01`。如果必须区分更高分辨率的蛋白等位型，输入完整等位基因名称，并确认参考表中存在该记录。

## 5. MHC-II 输入规则

MHC-II 最稳妥的输入方式是在同一个 `MHC` 单元格中同时提供 alpha 和 beta 等位基因。推荐使用 `/` 或 `+` 分隔：

```text
HLA-DRA*01:01/HLA-DRB1*04:01
HLA-DPA1*01:03+HLA-DPB1*04:01
HLA-DQA1*01:02/HLA-DQB1*06:02
```

alpha/beta 的先后顺序不影响输出；脚本会根据参考表中的链类型识别，并固定输出：

- `C`：MHC-II alpha；
- `D`：MHC-II beta。

默认的 `doc/mhcii_full_ipd/mhcii_model_input.csv` 已按成熟 MHC-II 结构域重新清理。
其中 `full_aa` 是去除信号肽、连接肽、跨膜区和胞质尾的可溶性结构预测用核心链，
不是完整蛋白前体。该表只保留蛋白序列无歧义的名称键；同一二字段名称对应多条
蛋白序列时不会静默选择第一条，而是从模型输入表排除。完整 QC 规则、来源记录和
冲突清单见 [`doc/mhcii_full_ipd/README.md`](doc/mhcii_full_ipd/README.md)。

配对时会检查：

- 必须恰好包含一条 alpha 和一条 beta；
- 两条链必须属于相同家族，即同为 `DR`、`DP` 或 `DQ`；
- 两条链必须来自相同物种。

### 单独输入 HLA-DRB

人源 HLA-DR 的 alpha 链多态性较低。只提供一个人源 `HLA-DRB` 时，脚本默认搭配：

```text
HLA-DRA*01:01
```

例如以下输入可直接转换：

```text
HLA-DRB1*04:01
```

可以修改默认 DRA：

```bash
python convert_tcr_pmhc_full_sequences.py input.csv \
  -o output.csv \
  --default-dr-alpha HLA-DRA*01:02
```

也可以禁止自动补充 DRA，强制每行显式提供 MHC-II 配对：

```bash
python convert_tcr_pmhc_full_sequences.py input.csv \
  -o output.csv \
  --default-dr-alpha ""
```

单独的 `DPA/DPB/DQA/DQB` 不能唯一决定另一条链，因此会被拒绝，必须显式输入 alpha/beta 等位基因对。

## 6. TCR 拼接逻辑

脚本不会使用简单的：

```text
V + CDR3 + J + C
```

因为 CDR3 的两端本身包含 V 和 J 来源的氨基酸，直接相加会重复残基，并且无法正确处理 V/J 末端切除。

脚本采用 trimming-aware 重建：

1. 在 V 基因末端窗口内，寻找与 CDR3 前缀相同的最长片段。
2. 在 J 基因起始窗口内，寻找与 CDR3 后缀相同的最长片段。
3. 根据匹配坐标裁剪 V、CDR3 和 J，使 CDR3 只出现一次。
4. 在 J 和常数区之间检测 suffix-prefix 重叠，再追加常数区。

核心重建形式为：

```text
V[:V_match_end] + CDR3[V_overlap_length:] + J[J_match_end:] + constant_region
```

常数区选择规则：

- TCR alpha：`TRAC`；
- `TRBJ1-*`：`TRBC1`；
- `TRBJ2-*`：`TRBC2`。

默认边界参数为：

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--min-overlap` | `3` | V/J 与 CDR3 至少需要匹配的氨基酸数 |
| `--max-overlap` | `30` | 搜索的最大匹配长度 |
| `--v-terminal-window` | `8` | V 匹配末端距离 V 参考序列末端的最大距离 |
| `--j-start-window` | `10` | J 匹配起点距离 J 参考序列起始端的最大距离 |

一般不建议为了提高通过率而随意放宽这些参数。匹配失败通常提示 V/J 注释、CDR3 定义、物种或参考序列之间存在不一致，应先检查输入数据。

## 7. 输出示例

输出 CSV 的表头固定为：

```csv
A,B,C,D,E
```

其中每一行与一条成功转换的输入记录对应：

```csv
<TCR_alpha_sequence>,<TCR_beta_sequence>,<MHC_alpha_sequence>,<B2M_or_MHC_beta_sequence>,<peptide_sequence>
```

在默认严格模式下，只有全部输入行均成功时才会写出结果，因此输入和输出行可以逐行对应。

## 8. 错误处理

### 严格模式：默认推荐

默认参数为：

```text
--on-error fail
```

只要任意一行出现问题，脚本就会停止，并打印 CSV 行号和错误原因；不会写出部分结果。适合正式结构预测和数据集构建。

```bash
python convert_tcr_pmhc_full_sequences.py input.csv -o output.csv
```

### 跳过错误行

需要先转换所有有效行时，可使用：

```bash
python convert_tcr_pmhc_full_sequences.py input.csv \
  -o output.csv \
  --on-error skip
```

脚本将生成：

```text
output.csv
output_rejected.csv
```

拒绝文件保留原始输入列，并增加：

| 列 | 含义 |
| --- | --- |
| `csv_row` | 原 CSV 中的实际行号，表头为第 1 行 |
| `conversion_error` | 拒绝原因 |

也可以指定拒绝文件路径：

```bash
python convert_tcr_pmhc_full_sequences.py input.csv \
  -o output.csv \
  --on-error skip \
  --rejected-output rejected_rows.csv
```

注意：使用 `skip` 后，输出行号不再与原输入行号一一对应；需要通过拒绝文件追踪被移除的行。

## 9. 参考序列与“全长”的含义

脚本重建的是“基于所选参考序列列的完整拼接链”。默认使用项目当前面向结构建模的序列列：

| 对象 | 默认参考列 |
| --- | --- |
| TCR V/J/C | `amino_acid` |
| MHC-I | `full_aa` |
| MHC-II | `full_aa` |

这些列是项目当前的 model-ready 输入，不应自动理解为一定包含信号肽、跨膜区和胞内尾的原始生物合成全长蛋白。对于 AlphaFold/tFold 的可溶性 TCR-pMHC 复合物预测，通常建议保留默认设置。

如确实需要参考表中的 original 序列，可以显式指定：

```bash
python convert_tcr_pmhc_full_sequences.py input.csv \
  -o output_original.csv \
  --gene-seq-column amino_acid_original \
  --mhci-seq-column full_aa_original \
  --mhcii-seq-column full_aa_original
```

切换为 original 序列前，应确认是否需要信号肽、跨膜区或胞内区。它们可能显著增加复合物预测长度，也可能影响 AlphaFold/tFold 的链构象和计算资源占用。

## 10. 自定义参考表

可以覆盖默认参考文件：

```bash
python convert_tcr_pmhc_full_sequences.py input.csv \
  -o output.csv \
  --gene-aa-path path/to/tcr_gene_aa.csv \
  --mhci-path path/to/mhci_sequences.csv \
  --mhcii-path path/to/mhcii_sequences.csv
```

参考表至少需要以下字段：

### TCR 表

```text
gene_name
<--gene-seq-column 指定的序列列>
```

### MHC-I 表

```text
MHC
<--mhci-seq-column 指定的重链序列列>
b2m_aa
```

### MHC-II 表

```text
name
chain
family
species
<--mhcii-seq-column 指定的序列列>
```

如果 MHC-II 表包含 `model_ready`，脚本只加载值为 `true/1/yes/y` 的记录。

## 11. 常见错误

### `CDR3A/CDR3B must start with the conserved C residue`

输入可能不是 IMGT 风格的完整 CDR3 氨基酸序列，或序列前端被截断。不要简单手动补 `C`；应先核对原始 V(D)J 注释。

### `CDR3 prefix has no trimming-aware overlap with TRAV/TRBV`

CDR3 前缀无法在指定 V 基因末端附近找到。常见原因包括：

- V 基因注释错误；
- CDR3 与 V 基因来自不同记录；
- 物种与参考表不一致；
- CDR3 定义不包含保守 C；
- 参考序列版本不同。

### `CDR3 suffix has no trimming-aware overlap with TRAJ/TRBJ`

CDR3 后缀与指定 J 基因不一致。应优先核对 J 基因和 CDR3 末端，而不是直接降低 `--min-overlap`。

### `MHC allele not found`

检查 MHC 命名、物种和参考表。脚本兼容 `HLA-A02:01` 与 `HLA-A*02:01`，但参考表中仍必须存在可匹配记录。

### `MHC-II pair must contain one alpha allele and one beta allele`

同一个 `MHC` 单元格中可能输入了两条 alpha 或两条 beta。请按 `DRA/DRB`、`DPA/DPB` 或 `DQA/DQB` 配对。

### `MHC-II family/species mismatch`

两条 MHC-II 链来自不同家族或不同物种。脚本会拒绝这种组合，避免生成不合理的异源配对。

## 12. 验证记录

开发时已用项目中的：

```text
data/epitope_tcr/pmhc_tcrab_fullseq_filtered.csv
```

进行全量回归验证。共比较 56,049 条记录，新脚本生成的 `A/B/C/D/E` 与原有 `tcra/tcrb/HLAa/b2m/peptide` 逐字符一致，五列差异数均为 0。

该结果验证了脚本与项目当前 MHC-I 全链转换主逻辑的一致性；新来源数据仍应检查拒绝记录和参考序列版本。
