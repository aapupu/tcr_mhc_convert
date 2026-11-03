# MHCtransform

`MHCtransform` is part of the [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) monorepo.

`MHCtransform` 是一个可审计、离线运行的人/鼠 MHC 转换项目，覆盖：

- MHC-I 重链与物种匹配的 B2M；
- MHC-II α/β 双链；
- IPD/UniProt deposited protein；其中通过官方 Full 状态或明确 UniProt 完整性检查的记录标为 `precursor_full`；
- 适合 AF3 复合物输入的成熟可溶性胞外链 `af3-soluble-extracellular-v1`；
- 两套严格区分的 34-aa pseudosequence；
- allele/name → sequence 与 sequence → **全部精确候选**的双向查询。

它不会把 pseudosequence 反向“猜成”一个唯一等位基因，也不会把去信号肽、去跨膜区的可溶性链称为全长。

## 表示定义

| 表示/`definition_id` | 含义 | 主要用途 |
|---|---|---|
| `precursor_full` | 来源记录中的完整前体：signal + extracellular + TM + cytoplasmic | 序列核对、表达构建设计的起点；不是默认 AF3 链 |
| `af3-soluble-extracellular-v1` | 去 signal/TM/cytoplasmic，保留全部胞外折叠域和膜近端连接肽 | AF3 pMHC–TCR 多链结构输入 |
| `netmhcpan-34-v1` | MHC-I 成熟重链参考对齐上的 34 个 peptide-contact 位点；不含 B2M | pan-allele 模型 descriptor，不是结构链 |
| `netmhciipan-34-v1` | MHC-II α 链 15 位 + β 链 19 位，按参考对齐拼接；允许 `X` | pan-allele 模型 descriptor，不是结构链 |

位点坐标和来源哈希保存在 [`metadata.json`](src/mhctransform/data/metadata.json)。MHC-II 不能直接用统一的 precursor 下标切片：信号肽长度和对齐插缺不同，因此运行时使用带版本的已发布 pseudosequence artifact，而不是从全长链临时截取。

## 当前本地参考库

| 内容 | human I | human II | mouse I | mouse II |
|---|---:|---:|---:|---:|
| deposited protein records | 20,066 | 9,419 | 12 | 17 |
| 默认可用的 complete precursor records | 11,283 | 2,525 | 9 | 14 |
| AF3-ready extracellular-chain records | 11,004 | 2,456 | 9 | 14 |
| 34-aa records | 2,953 | 5,620 | 7 | 8 |

人 HLA 来自固定的 IPD-IMGT/HLA 3.66.0 / commit `b2055ea`，并与同版 `Allele_status.txt` 的官方 `Full/Partial` 状态联合判定；鼠 MHC 和人/鼠 B2M 来自显式 UniProt accession。AF3 边界由 reviewed UniProt topology 投影到同版 IPD 蛋白比对，且逐条验证为前体的连续子序列。34-aa 科学定义来自 NetMHCpan/NetMHCIIpan 文献，本地分发数据采用 PUFFIN 的公开二级镜像。完整来源、accession、内容 SHA-256 和生成日期均在 metadata 中。

deposited records 包含不完整条目；默认查询只返回官方状态为 `Full` 且没有未知残基的人链，或通过 accession/fragment/长度检查的鼠链。`--include-partial` 才会显式纳入片段。当前鼠库是常用、可溯源等位型集合，不声称穷尽全部 H-2 等位基因。

## 使用

无需第三方运行依赖：

```bash
cd MHCtransform
python3 -m pip install -e .
```

### 等位基因到前体全长

```bash
mhctransform allele-to-full 'HLA-A*02:01:01:01' --class I
mhctransform allele-to-full 'H-2-Kb' --class I
mhctransform allele-to-full 'HLA-DRB1*15:01' --class II
```

两字段名称仍可能对应多个完整蛋白。例：`HLA-C*03:03` 在本版数据中有 366 aa 与 372 aa 两个官方 Full 候选，会标记 `exact_ambiguous`；高字段名称可进一步消歧。`HLA-A*02:01` 的旧 364-aa 条目按官方状态为 Partial，默认只返回 365-aa Full 候选。

不存在的高字段名称默认不会静默降到两字段。例如 `HLA-A*02:01:999:999` 返回 `not_found`。只有调用者明确给出 `--allow-two-field-projection` 时才投影，并用 `projected_unique/projected_ambiguous` 标记，绝不伪装成 exact。官方已知但只有 Partial protein 的精确名称即使开启该开关也不投影。HLA `N/L/S/Q/C/A` 表达后缀及其完整官方字段被保留；裸 `HLA-A*03:03` 不会误命中 `HLA-A*03:03N`，也不会把官方 `HLA-A*02:01:01:02L` 人工改成不存在的 `HLA-A*02:01L`。

### 完整 MHC 分子到两条链

```bash
# MHC-I：重链 + 物种匹配 B2M
mhctransform molecule-to-chains 'H-2-Kb' --class I

# MHC-II：显式 α/β 配对
mhctransform molecule-to-chains \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# 鼠分子别名
mhctransform molecule-to-chains 'H-2-IEk' --class II
```

仅给 human DRB 时，程序可按常见约定补 `HLA-DRA*01:01`，但结果的 `notes` 会明确披露这一假设。DQ/DP 必须给出 α/β 配对；程序不仅拒绝无约束笛卡尔积，也会拒绝 DQα/DPβ、DPα/DQβ 等跨家族杂配。

### AF3 胞外链与两链配对

```bash
# 单链：allele -> AF3 chain
mhctransform allele-to-af3 'HLA-A*02:01' --class I

# 反向：AF3 chain -> 全部精确候选
mhctransform af3-to-allele @af3_chain.txt --class I --species human

# 完整 MHC-I 分子：heavy 284 aa + human B2M 99 aa，两条独立链
mhctransform molecule-to-af3 'HLA-A*02:01' --class I

# 完整 MHC-II 分子：α、β 两条独立链
mhctransform molecule-to-af3 \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II
```

`molecule-to-af3` 不会把两条 MHC 链拼接，也不添加 linker/tag。MHC-I 自动加入同物种 B2M；MHC-II 要求同家族 α/β。默认 AF3 registry 排除 partial 记录以及带 `N/L/S/Q/C/A` 表达后缀的人 HLA，并在能定位原因时通过 `notes` 返回排除理由。

这里的裁剪是项目的可审计工程 profile，不是 Google 官方指定的 pMHC–TCR 坐标。AF3 官方只规定把各分子作为独立 entity 并提供 amino-acid sequence；34-aa pseudo 是非连续 descriptor，严禁作为 AF3 蛋白链。

### 34 aa 双向转换

```bash
mhctransform allele-to-pseudo 'HLA-A*02:01' \
  --class I --definition netmhcpan-34-v1

mhctransform pseudo-to-allele \
  YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY \
  --class I --definition netmhcpan-34-v1

mhctransform allele-to-pseudo 'HLA-DRB1*15:01' \
  --class II --definition netmhciipan-34-v1

mhctransform pseudo-to-allele \
  YSYFLASGGQVVHVLYFGYTYHDIRTETVHGPHT \
  --class II --definition netmhciipan-34-v1
```

反向转换只做 exact lookup。一个 34mer 常由许多等位基因共享，因此结果通常是 `exact_ambiguous`；未知 34mer 返回 `not_found`，不会偷偷做最近邻预测。

### 全长序列反查

```bash
mhctransform full-to-allele @sequence.txt --class I --species human
```

输出中的一个 sequence candidate 可能把多个完全相同的 deposited HLA 名称压缩在 `aliases` 内；此时状态仍会标为 `exact_ambiguous`。

### Python API

```python
from mhctransform import MHCRegistry
from mhctransform.registry import AF3_DEFINITION, PSEUDO_I_DEFINITION

registry = MHCRegistry()

forward = registry.allele_to_pseudo(
    "HLA-A*02:01",
    mhc_class="I",
    definition_id=PSEUDO_I_DEFINITION,
)

reverse = registry.pseudo_to_alleles(
    forward.candidates[0].pseudo_sequence,
    mhc_class="I",
    definition_id=PSEUDO_I_DEFINITION,
)
print(reverse.status, reverse.candidate_count)

af3_molecule = registry.molecule_to_af3_chains(
    "HLA-A*02:01",
    mhc_class="I",
    definition_id=AF3_DEFINITION,
)
for chain in af3_molecule.candidates[0].chains:
    print(chain.chain_role, len(chain.sequence), chain.sequence)
```

## 金标准结果

| 查询 | precursor | AF3 profile（前体 1-based） | 34 aa |
|---|---:|---|---|
| `HLA-A*02:01` | 365 aa | 25–308，284 aa，`GSHSMRY...QPTIPI` | `YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY` |
| `B2M-human` | 119 aa | 21–119，99 aa，`IQRTPK...KWDRDM` | 不适用 |
| `H-2-Kb` | 369 aa | 22–305，284 aa，`GPHSLR...STVSNM` | `YVEYYREKAGNSFVDTLYIVSQYYTWAELAYTWY` |
| `HLA-DRA*01:01` | 254 aa | 26–216，191 aa | 与下一行组合 |
| `HLA-DRB1*15:01` | 266 aa | 30–227，198 aa | `QEFFIASGAAVDAIMWPRFDYFDIQAATYHVVFT` |
| `H-2-IAb` α/β | 256/265 aa | 24–218（195 aa）/ 28–226（199 aa） | `YSYFLASGGQVVHVLYFGYTYHDIRTETVHGPHT` |

HLA-A*02:01 的 AF3 链最后是 `...QPTIPI`；下一位 precursor 309 的 `V` 已进入 UniProt 标注的 TM 309–332，因此不能多保留。旧实现的 274-aa HLA-I 裁剪则少了 299–308 的膜近端 connector，也不作为默认 AF3 输出。

## 重新生成数据

生成脚本会复用本地缓存；首次运行需要访问固定版 IPD-IMGT/HLA、PUFFIN 与 UniProt：

```bash
python3 tools/build_reference.py
python3 tools/build_reference.py --offline
```

脚本不会从旧版 `mhcii_model_input.csv` 继承裁剪结果。人 HLA full registry 用同版 FASTA 和 `Allele_status.txt` 联合构建；AF3 链用固定版 IPD protein alignment 把 exact allele 映射到 reviewed UniProt extracellular topology，并要求结果可唯一定位回 parent precursor；pseudo artifact 再用官方 2009→现行命名表和 deleted/renamed 表规范旧名。鼠链按显式 accession 获取；可变的 UniProt REST 批次另外用预期 SHA-256 固定，内容变化会 fail closed，不能静默改变边界。所有生成文件均稳定排序并写入来源哈希。二级 pseudo 镜像中的 5 个源冲突没有静默保留或任选首条，而是按可审计规则写入 `quarantined_pseudosequences.csv`。

## 测试

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

28 项测试覆盖人/鼠、MHC-I/II、full/pseudo/AF3 正反向、12 个 AF3 金标准边界、13,483 条 AF3 链的 parent-subsequence invariant、MHC-II 家族配对、B2M、`X` 位点、官方 Partial fail-closed、表达后缀隔离、多对一歧义、pre-2010 HLA 名称、三位 DPB1 字段、compact-name 碰撞和显式 projection 状态。

## 重要边界

- `precursor_full` 表示来源确认覆盖完整 coding protein 的记录，通常包含 signal peptide、extracellular domains、TM 和 cytoplasmic tail；真实 null/截短等位基因也可能是完整记录。它不是 AF3 默认输入。
- AF3 的结构输入可直接使用 `af3-soluble-extracellular-v1`；不要把 34-aa pseudosequence 当作可折叠链。
- P/G group 只表达特定结构域的一致性，不能当作整条前体蛋白完全一致。
- 当前 pseudosequence artifact 与最新 IPD 全长快照覆盖范围不同；没有 exact pseudo 记录时返回 `not_found`。
- 源数据的许可与使用条款不因本项目代码许可而改变，见 [`NOTICE.md`](NOTICE.md)。

详细审核、TCR 工具对比和 AF3 五链建议见 [`AUDIT_REPORT.md`](AUDIT_REPORT.md)。

