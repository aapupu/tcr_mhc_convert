# TCRtransform

[English](README.md) | **简体中文**

`TCRtransform` 是 [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) 单一仓库中的独立 Python 包。

`TCRtransform` 是一个可独立安装的 Python 工具，用于人和小鼠 αβ TCR 的：

- 从 `V gene + CDR3 + J gene + constant` 重建结构预测用完整胞外氨基酸链，可用于AF3输入；
- 根据 V 基因注释 CDR1、CDR2 和 CDR2.5；
- 定义 CDR3 内 V 支持前缀、J 支持后缀和中间的非 V/J 可解释区间；
- 从完整 TCRα/TCRβ 氨基酸链反向推断 V/J/C、CDR1/2/2.5/3 及 CDR3 边界；
- 对 CSV 进行可重复、原子化的批量处理。

仅依赖 Python 标准库，支持 Python 3.9 及以上版本。

## 生物学定义

### CDR1、CDR2 和 CDR2.5

CDR1、CDR2、CDR2.5 来自打包的 IMGT 官方 V/J 快照，使用与 tcrdist 兼容的对齐列（包括小鼠 TRA 的额外列）。完整链、CDR 和 junction 元数据使用同一次解析得到的等位基因。CDR2.5 是 tcrdist3 所称的 pMHC-facing loop/HV4 区域，并不是额外的一条经典 IMGT CDR。

### CDR3

本项目接受和输出包含两端锚点的 CDR3 氨基酸序列：从 V 侧保守 `C` 开始，到 J 侧注释锚点结束。J 锚点通常为 `F/W`，也支持参考表中明确记录的非经典锚点。

CDR3 被描述为：

```text
[V-supported prefix][protein-level non-V/J interval][J-supported suffix]
```

`CDR3*_non_vj` 只表示未被所选 V/J 蛋白参考覆盖的氨基酸区间。它不能单独证明真实的核苷酸 N/P 插入来源，原因包括密码子分相、趋同氨基酸和 β 链 D 基因贡献。若需要核苷酸层面的精确归因，应使用完整 V(D)J 核苷酸序列及专门的重排注释工具。

α 链和 β 链都允许出现 `non_vj_length = 0`。工具不会为了给优化模型留下可编辑位点而强制修改 germline 支持残基。

## 拼接原则

TCR 链不是简单执行 `V + CDR3 + J + C`。工具会：

1. 在 V 末端窗口中寻找 CDR3 前缀的最长 trimming-aware 匹配；
2. 在 J 起始窗口中寻找 CDR3 后缀的最长 trimming-aware 匹配；
3. 裁掉重复部分后拼接，保证 CDR3 只出现一次；
4. β 链根据 `TRBJ1-*`/`TRBJ2-*` 选择 `TRBC1`/`TRBC2`；
5. α 链根据参考元数据补回 TRAJ–TRAC 分相密码子产生的桥接氨基酸；
6. 追加结构预测用常数区胞外结构域。

默认 `min_overlap=3` 是严格的重建质控阈值。若用户显式设置为 1 或 2，短匹配必须同时满足：

- V 侧从保守 `C` 开始；
- J 侧终止于该 J 的注释锚点；
- 匹配位于参考记录的预期 V/J 边界；
- 在允许窗口内位置唯一。

默认阈值拒绝 1–2 aa 匹配只代表证据不足，不代表该 TCR 在生物学上错误。

## 安装

```bash
cd TCRtransform
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

安装后可使用 `tcrtransform`，也可以直接运行 `python -m tcrtransform`。

## CSV 命令

### 1. 完整拼接并注释全部区域

输入至少包含 `TRAV, TRAJ, CDR3A, TRBV, TRBJ, CDR3B`：

```bash
tcrtransform build examples/tcr_input.csv \
  --species-column species \
  -o examples/tcr_built.csv
```

输出保留原列，并增加 `TCRa`、`TCRb`、CDR1/2/2.5、CDR3 边界、固定掩码、参考解析方式和状态列。

### 2. 只注释 CDR 区域

```bash
tcrtransform annotate examples/tcr_input.csv \
  --species-column species \
  -o examples/tcr_annotated.csv
```

### 3. 从完整链反向解析

输入包含 `TCRa/TCRb`，也接受 `A/B` 等常见别名：

```bash
tcrtransform parse full_tcr.csv \
  --species human \
  -o full_tcr_parsed.csv
```

若每行物种不同，可改用 `--species-column species`。默认遇到任意错误时不写最终输出；使用 `--on-error keep` 可保留失败行并在 `tcrtransform_status` 中记录原因。

输出文件已存在时需显式添加 `--overwrite`。输入中已有目标注释列时，默认拒绝覆盖；确认需要重新计算时添加 `--overwrite-columns`。

## 关键输出列

以 α 链为例，β 链使用相同字段并以 `b` 结尾：

| 列 | 含义 |
|---|---|
| `CDR1a`, `CDR2a`, `CDR2.5a`, `CDR3a` | CDR 氨基酸序列 |
| `CDR3a_v_supported` | 与 V 参考边界一致的 CDR3 前缀 |
| `CDR3a_non_vj` | 未被已接受 V/J 匹配覆盖的连续区间 |
| `CDR3a_j_supported` | 与 J 参考边界一致的 CDR3 后缀 |
| `CDR3a_keep_prefix`, `CDR3a_keep_suffix` | 优化时建议固定的两端长度 |
| `CDR3a_fixed_mask` | 与 CDR3 等长，`1` 为 V/J 支持位点，`0` 为未覆盖位点 |
| `CDR3a_non_vj_start/end` | 0-based、左闭右开的候选可编辑区间 |
| `CDR3a_safe_for_optimization` | V、J 两侧是否都成功解析 |
| `CDR3a_boundary_status` | `ok`、`fully_vj_supported`、短匹配低于阈值或未解析等状态 |
| `TCRa_reference_resolution` | V/J/C 重建使用精确参考还是 gene-level fallback |

下游优化应至少要求 `safe_for_optimization=True`，并直接使用 `fixed_mask`。不要通过强制 `min_editable > 0` 覆盖已经被 V/J 支持的位点。

## Python API

```python
from tcrtransform import (
    annotate_v_gene,
    build_paired_tcr,
    define_cdr3_boundary,
    parse_chain,
)

loops = annotate_v_gene("human", "B", "TRBV20-1")
print(loops.cdr1, loops.cdr2, loops.cdr2_5)

boundary = define_cdr3_boundary(
    species="human",
    chain="B",
    v_gene="TRBV30",
    j_gene="TRBJ1-1",
    cdr3="CAWGRGARAFIDTEAFF",
)
print(boundary.non_vj_region, boundary.fixed_mask)

pair = build_paired_tcr(
    species="mouse",
    trav="TRAV14D-1",
    traj="TRAJ33",
    cdr3a="CAASDNYQLIW",
    trbv="TRBV12-1",
    trbj="TRBJ2-7",
    cdr3b="CASSRANYEQYF",
)
print(pair.alpha.sequence)
print(pair.beta.sequence)

parsed = parse_chain(pair.alpha.sequence, "A", species="mouse")
print(parsed.v_genes, parsed.j_genes, parsed.cdr3)
```

## 参考数据与适用范围

- 支持物种：`Homo sapiens` 和 `Mus musculus`。
- 优先使用已收录的 V/J 等位基因（`*01`、`*02`、`*03` 及其他可用编号）。未提供编号时使用 gene-level 代表序列；编号未收录或功能性等位基因序列不完整时，回退到同一 gene-level 代表。基因也不存在则报错。`TCR*_reference_resolution` 按 V/J/C 分别标记 `exact`、`gene_level`、`gene_level_fallback`。C 区仍使用 gene-level 参考。
- Germline 的 P/ORF 标签或终止符不能单独用于拒绝实际重排的受体。有可靠边界的参考先进行拼接；在有证据支持的 V/J 接头中被移除的终止符允许通过并记录，拼接后仍保留终止符则报错并要求复核。绝不直接删除保留区域内的终止符。P/ORF 作为复核提示，不代表样本一定不 productive。
- `TCR*_reference_warnings` 记录 germline 标签与被移除的终止符。`TCR*_productivity_status` 为 `not_assessed`，因为氨基酸重建不能证明阅读框、剪接完整性和生物学 productive。CSV 默认 `--on-error keep`，保留失败行与原始 `productive` 列。剔除或重新生成训练/测试数据前，应结合实际重排序列、保留的终止符位置、核苷酸证据与原始 productive 标注。仍可显式使用 `--on-error fail`。
- 官方快照含人/鼠 V/J 共 870 条记录。历史 `usable` 标记表示最初的 684 条严格子集（人：263 V + 76 J；鼠：277 V + 68 J），不表示当前重排候选资格；运行时也考虑具有可靠边界的 P/ORF 和含终止符参考。gene-level 输入不能恢复未提供的个体多态性；忠实克隆需要实测序列或自定义参考。
- 默认常数区为结构预测用胞外结构域，不包含 leader、跨膜区和胞内尾。因此这里的“完整链”指本项目的结构模型输入链，不等同于完整表达载体 CDS。
- 反向解析保留蛋白序列无法区分的并列 V/J 基因，并用 `|` 连接，不会任意挑选一个。
- 不支持 γδ TCR。

参考方法与定义：

- [Stitchr, Nucleic Acids Research 2022](https://doi.org/10.1093/nar/gkac190)
- [IMGT FR/CDR definition](https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html)
- [tcrdist3 database conventions](https://github.com/kmayerb/tcrdist3/blob/master/docs/inputs.rst)

数据快照、commit 和校验值见 `src/tcrtransform/data/README.md` 与 `NOTICE.md`。

`python tools/build_imgt_reference.py` 从 IMGT 官方更新快照及派生的 gene/junction 表；加 `--from-bundled` 可离线重建。旧的 `build_mouse_reference.py` 和 `build_junction_reference.py` 仅代表历史 tcrdist 构建流程，不应用于覆盖已修正的默认参考。数据目录保留来源 URL、accession、功能性及下载哈希。

例如人 `TRAJ24*01` → `TTDSWGKFEFGAGTQVVVTP`，`TRAJ24*02` → `TTDSWGKLQFGAGTQVVVTP`，`TRAJ24*03` → `TTDSWGKFQFGAGTQVVVTP`。`TRAJ24` 与未收录的 `TRAJ24*99` 均使用 `*02` gene-level 代表。原错误序列 `TTDSWGKLQFEFGAGTQVVVTP` 已停用。小鼠 `TRAV14-1`/`TRAV14D-1`、`TRAV13D-1` 及 J 区修正见[修正报告](validation/IMGT_CORRECTIONS.md)。

## 测试

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
```

测试包含人/鼠 CDR 参考、短匹配防误判、完全由 V/J 支持且可编辑区长度为 0、TRAJ–TRAC 桥接、CSV 处理和正反向 round-trip。

## 对照验证资料

`validation/` 保留了本包核心实现从原项目抽取前进行的随机对照实验：人和小鼠各 100 条完整双链 TCR，与 tcrdist3 和 Stitchr 对照。摘要见 [`validation/REPORT.md`](validation/REPORT.md)，逐链结果、失败分类、抽样记录和机器可读汇总也一并保留，便于审计。

## 项目结构

```text
TCRtransform/
├── pyproject.toml
├── README.md
├── README.zh-CN.md
├── NOTICE.md
├── LICENSE
├── examples/
├── validation/
├── tools/
├── src/tcrtransform/
│   ├── cli.py
│   ├── junctions.py
│   ├── references.py
│   ├── regions.py
│   ├── reverse.py
│   ├── transform.py
│   └── data/
└── tests/
```

## 许可证与数据条款

本包代码采用 [MIT License](LICENSE)。捆绑及派生的参考数据仍保留上游署名和适用条款；MIT 许可证不会重新许可 IMGT/tcrdist 派生数据。
