# MHCtransform

[English](README.md) | **简体中文**

`MHCtransform` 是 [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) 单一仓库中的独立 Python 包。

`MHCtransform` 是一个可审计、离线运行的人/鼠 MHC 转换项目，覆盖：

- MHC-I 重链与物种匹配的 B2M；
- MHC-II α/β 双链；
- IPD/UniProt deposited protein；其中通过官方 Full 状态或明确 UniProt 完整性检查的记录标为 `precursor_full`；
- 适合 AF3 复合物输入的成熟可溶性胞外链 `af3-soluble-extracellular-v1`；
- 连续肽结合区：MHC-I 输出一条 alpha1+alpha2，MHC-II 输出独立的 alpha1 与 beta1 两条链；
- 与严格参考库隔离的 partial、partial 引导补全和 pseudo-only scaffold 后备层；
- 两套严格区分的 34-aa pseudosequence；
- 双向查询：严格序列返回**全部精确候选**；显式纳入的 fallback 构建只返回 `compatible_*` 候选。

它不会把 pseudosequence 反向“猜成”一个唯一等位基因，也不会把去信号肽、去跨膜区的可溶性链称为全长。

## 表示定义

| 表示/`definition_id` | 含义 | 主要用途 |
|---|---|---|
| `precursor_full` | 来源记录中的完整前体：signal + extracellular + TM + cytoplasmic | 序列核对、表达构建设计的起点；不是默认 AF3 链 |
| `af3-soluble-extracellular-v1` | 去 signal/TM/cytoplasmic，保留全部胞外折叠域和膜近端连接肽 | AF3 pMHC–TCR 多链结构输入 |
| `source-partial-v1` | 精确 deposited fragment，绝不重新标记为完整前体 | 证据恢复和片段检查 |
| `precursor-imputed-v1` | 合成完整长度 scaffold：保留对齐到的 deposited 残基及目标 pseudo 的全部非 `X` 残基，其余来自同物种、同位点模板 | 显式开启的构建设计；不是官方等位基因序列 |
| `af3-soluble-extracellular-imputed-v1` | 使用相同逐残基来源规则生成的 AF3 长度构建链 | 没有 exact 链时显式开启 AF3 初筛 |
| `mhci-alpha1-alpha2-v1` | 连续 MHC-I alpha1+alpha2 肽结合槽，通常为 182 aa；不含 alpha3 与 B2M | 结构域比较及 34 个接触位点的可审计提取/核验 |
| `mhcii-alpha1-v1` | 连续 MHC-II alpha1 肽结合域 | MHC-II PBR 的第一条链 |
| `mhcii-beta1-v1` | 连续 MHC-II beta1 肽结合域 | MHC-II PBR 的第二条链 |
| `netmhcpan-34-v1` | MHC-I 成熟重链参考对齐上的 34 个 peptide-contact 位点；不含 B2M | pan-allele 模型 descriptor，不是结构链 |
| `netmhciipan-34-v1` | MHC-II α 链 15 位 + β 链 19 位，按参考对齐拼接；允许 `X` | pan-allele 模型 descriptor，不是结构链 |

PBR 与 pseudo 相互配合但不能互换：PBR 是连续、可折叠的结构域并包含全部 pseudo 接触位点；34mer 是非连续模型 descriptor。PBR 适合比较与 pseudo 核验，但 pMHC–TCR 的 AF3 默认仍使用完整胞外链，以保留膜近端结构域及天然界面。

位点坐标和来源哈希保存在 [`metadata.json`](src/mhctransform/data/metadata.json)。MHC-II pseudo 使用经过校准的参考对齐列，不能把统一 precursor 下标直接套到所有链。后备策略逐级累积且默认关闭：`none` → `partial` → `imputed` → `pseudo-scaffold`。34mer 仍是多对多：15,478 个人 MHC-I 名称仅对应 4,279 个唯一 pseudo，8,948 个人 MHC-II 分子名称仅对应 1,695 个。因此 imputed 输出只是兼容 scaffold，不是真实全长恢复。

## 当前本地参考库

| 内容 | human I | human II | mouse I | mouse II |
|---|---:|---:|---:|---:|
| deposited protein records | 20,066 | 9,419 | 12 | 17 |
| 默认可用的 complete precursor records | 11,283 | 2,525 | 9 | 14 |
| 严格 AF3-ready extracellular-chain records | 11,004 | 2,456 | 9 | 15 |
| 严格 PBR chain records | 11,003 | 2,456 | 8 | 15 |
| 34-aa records | 15,478 | 8,948 | 10 | 8 |

人 HLA 来自固定的 IPD-IMGT/HLA 3.66.0 / commit `b2055ea`，并结合同版官方 `Full/Partial` 状态；鼠 MHC 和人/鼠 B2M 来自显式 UniProt accession。AF3 与 PBR 边界由 reviewed UniProt topology/domain 标注投影到同版 IPD 蛋白比对，严格输出逐条验证为 parent precursor 的连续子序列。

34-aa registry 以 NetMHCpan 4.2f 与 NetMHCIIpan 4.3k 的公开支持名称快照为标准。人/鼠范围分别含 10,395 个 MHC-I 名称和 11,344 个 MHC-II 分子；本版分别为其中 8,743 与 8,956 个提供经审查的 34-aa 伪序列。出现在支持清单中并不等于已知序列：若 deposited source 与经审查 legacy 行都不能覆盖所需接触序列，构建器会报告不可用，不会补猜 pseudo。

捆绑的 24,444 行 pseudo（MHC-I 15,488；MHC-II 8,956）明确区分证据：`derived_exact_pbr` 13,560、`derived_contact_complete_full` 89、`derived_contact_complete_partial` 7,223、`legacy_secondary` 3,572；没有捆绑任何 `official_reference` 行。可选的本地 DTU training pseudo 仅作为 validation oracle，不会把 oracle-only 行复制进 registry。Oracle 精确核验为 MHC-I 215/216（只缺 `H-2-Kq`）和 MHC-II 5,565/5,631。4 条含内部 indel 的 MHC-I 记录通过 oracle 在独立 alignment contact-column 候选间进行选择后安全恢复；其 sequence evidence 仍是 deposited protein，不是复制 oracle。partial 只有在全部接触位点无歧义映射且没有 `X`、gap 或插补时才能生成 pseudo。

现行命名处理也显式记录：历史 `HLA-DQA1*01:07` 与 `HLA-DRB1*15:13` 输入分别解析到当前 `HLA-DQA1*01:07Q` 与 `HLA-DRB1*15:13Q`。`HLA-DRB1*13:10` 保留的 `legacy_secondary` 行是 PUFFIN/NetMHC-compatible `ESSFDSYSIDKATYHVGFT`（`official_training_exact`）；它带有与 current-IPD-derived `ESSFDYFDIDKATYHVVFT` 冲突的质量标志，而不是静默覆盖任一观察。Component registry 含 4,386 条链 pseudo，另有 280 条拒绝/冲突记录隔离。

`pseudo_components.csv` 是只包含已解析、可独立证明链组分的审计表，并不保证每条 MHC-II pair 都能完整拆分。具体而言，64 条保留的 `legacy_secondary` pair（4 个有争议的 DPB1 beta 组分 × 16 个 DPA1 alpha partner）因 legacy、current-IPD-derived 与本地 oracle 证据冲突，故意不把未解析 beta 写入 component 表。Pair 行连同来源仍保留，但争议 beta 序列不会被提升为独立组分。

deposited records 包含不完整条目。默认查询只返回官方 `Full` 且无未知残基的人链，或通过 accession/fragment/长度检查的鼠链；片段与合成构建必须显式开启。`partial` 不会补尾，`imputed` 要求 deposited 引导，`pseudo-scaffold` 是最低证据级别。带 `N/L/S/Q/C/A` 表达后缀的人 HLA 默认仍禁止进入 AF3。当前鼠库是常用、可追溯 H-2 集合，不声称穷尽全部 IPD-MHC。

重建后的 fallback registry 含 18,079 个 component targets、13,217 条 full/fragment、6,626 条 AF3 和对应的 6,626 条 PBR。它们是构建/来源片段数量，不代表新增了同等数量的已验证生物学全长序列。

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

# 仅返回精确 deposited fragment；representation=source-partial-v1
mhctransform allele-to-full 'HLA-B*08:03' --class I --fallback partial

# 用 deposited fragment 约束的完整长度合成构建
mhctransform allele-to-full 'HLA-B*08:03' --class I --fallback imputed
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

### 肽结合区链

```bash
# MHC-I：一条连续 alpha1+alpha2 重链区，通常 182 aa
mhctransform allele-to-pbr 'HLA-A*02:01' --class I
mhctransform molecule-to-pbr 'HLA-A*02:01' --class I

# MHC-II：独立输出 alpha1 与 beta1 两条链
mhctransform molecule-to-pbr \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# 反查返回所有精确链候选
mhctransform pbr-to-allele @pbr_chain.txt \
  --class I --definition mhci-alpha1-alpha2-v1
```

MHC-I 的 `molecule-to-pbr` 有意省略 alpha3 与 B2M；MHC-II 返回两个 sequence object，绝不人工拼接。常见 reviewed domain 长度为 MHC-I 182 aa、MHC-II alpha1 84–96 aa、beta1 90–95 aa；实际边界按位点和等位基因投影。后备 PBR 与 AF3 后备一样保留 observed/imputed spans 和证据等级。

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

### Partial 与 imputed AF3 后备

```bash
# 不输出长序列，一次查看各表示层是否可用
mhctransform coverage 'HLA-B*08:03' --class I --species human

# 优先 partial 引导构建；若只有 pseudo-only scaffold 则仍失败
mhctransform molecule-to-af3 'HLA-B*08:03' \
  --class I --fallback imputed

# 没有 deposited protein 时允许最低置信度 pseudo-only scaffold
mhctransform molecule-to-af3 'HLA-A*01:04' \
  --class I --fallback pseudo-scaffold

# MHC-II 仍输出两条独立链
mhctransform molecule-to-af3 'H-2-IAq' \
  --class II --fallback pseudo-scaffold
```

每个后备候选都会报告 `evidence_level`、deposited-observed/imputed 区间与数量、来源记录、模板记录/等位基因、pseudo 距离、未确定的 `X` 位点及质量标志。若 deposited fragment 已覆盖完整胞外目标，则提升为 `source_exact_ectodomain`；更短的片段绝不直接作为 AF3 链，必须进入 `imputed` 模式。未直接观察的位置包含模板来源与 pseudo 约束残基，pseudo 来源另行保留。`H-2-IAs` 不再使用推断：UniProt P14437 提供精确的 1–195 α 链胞外区，可与精确 β 链配对；`H-2-IAq` 与 `H-2-IAu` 的 α 链仍需 imputed。反向查询默认排除后备构建；只有显式给出 `--include-fallback` 才纳入，并只表示构建兼容性，不证明真实生物学等位基因。

审核集中稀有名称的代表性结果：

| 名称 | deposited partial | AF3 输出 | 证据/模板 |
|---|---:|---:|---|
| `HLA-B*08:03` | 181 aa | 285-aa heavy + B2M | partial-guided；`HLA-B*08:02` |
| `HLA-B*45:06` | 181 aa | 285-aa heavy + B2M | partial-guided；`HLA-B*45:01` |
| `HLA-A*02:50` | 181 aa | 284-aa heavy + B2M | partial-guided；`HLA-A*02:02` |
| `HLA-A*03:19` | 181 aa | 284-aa heavy + B2M | partial-guided；`HLA-A*03:01` |
| `H-2-Lq` | 206 aa | 285-aa heavy + mouse B2M | partial-guided；`H-2-Ld` |
| `H-2-IAs` | α 为 233-aa mature membrane record | α 195 aa + β 197 aa | 两条均 strict exact；α 来自 P14437 |
| `H-2-IAq` | α fragment | α 195 aa + β 200 aa | α 由 `H2-Aa*k` 补全；β exact |
| `H-2-IAu` | α fragment | α 195 aa + β 197 aa | α 由 `H2-Aa*k` 补全；β exact |

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

反向转换只做 exact lookup。一个 34mer 常由许多等位基因共享，因此结果通常是 `exact_ambiguous`；未知 34mer 返回 `not_found`，不会偷偷做最近邻预测。MHC-II 34mer 表示一个 alpha/beta 分子：先 15 个 alpha 位点，再 19 个 beta 位点。支持空间中的 DQ/DP 组合是模型输入定义，不证明常见 haplotype 或体内配对。

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
    "HLA-B*08:03",
    mhc_class="I",
    definition_id=AF3_DEFINITION,
    fallback="imputed",
)
for chain in af3_molecule.candidates[0].chains:
    print(chain.chain_role, len(chain.sequence), chain.sequence)
```

## 金标准结果

| 查询 | precursor | AF3 profile（前体 1-based） | PBR | 34 aa |
|---|---:|---|---|---|
| `HLA-A*02:01` | 365 aa | 25–308，284 aa，`GSHSMRY...QPTIPI` | alpha1+alpha2，182 aa | `YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY` |
| `B2M-human` | 119 aa | 21–119，99 aa，`IQRTPK...KWDRDM` | 不适用 | 不适用 |
| `H-2-Kb` | 369 aa | 22–305，284 aa，`GPHSLR...STVSNM` | alpha1+alpha2，182 aa | `YVEYYREKAGNSFVDTLYIVSQYYTWAELAYTWY` |
| `HLA-DRA*01:01` | 254 aa | 26–216，191 aa | alpha1，84 aa | 与下一行组合 |
| `HLA-DRB1*15:01` | 266 aa | 30–227，198 aa | beta1，95 aa | `QEFFIASGAAVDAIMWPRFDYFDIQAATYHVVFT` |
| `H-2-IAb` α/β | 256/265 aa | 24–218（195 aa）/ 28–226（199 aa） | 独立 alpha1/beta1 | `YSYFLASGGQVVHVLYFGYTYHDIRTETVHGPHT` |

HLA-A*02:01 的 AF3 链最后是 `...QPTIPI`；下一位 precursor 309 的 `V` 已进入 UniProt 标注的 TM 309–332，因此不能多保留。旧实现的 274-aa HLA-I 裁剪则少了 299–308 的膜近端 connector，也不作为默认 AF3 输出。

## 重新生成数据

核心生成脚本会复用本地缓存；首次在线运行需要访问固定版 IPD-IMGT/HLA 与 UniProt。没有本地 DTU oracle 文件时仍可重建独立 strict/PBR 与 secondary 资产，但不能复现经过 oracle gate 的发布 registry：

```bash
python3 tools/build_reference.py
python3 tools/build_reference.py --offline

# 仅从已有严格资产重建确定性的 fallback 数据
# 同时重建 fallback_pbr_sequences.csv
python3 tools/build_fallback_reference.py

# 需要时直接重建 strict 与 fallback PBR
python3 tools/build_pbr_reference.py
```

DTU 的公开支持名称快照以固定 SHA-256 保存在 `tools/reference/`。要复现发布数据，可把 training pseudo 下载到本地缓存，并且只作为 validation oracle 传给集成 builder；原始归档、解出的 oracle 文件和 oracle-only 行都不会提交：

```bash
python3 tools/fetch_netmhc_training_pseudos.py --cache-dir /path/to/netmhc-cache

python3 tools/build_reference.py \
  --netmhcpan-training-pseudo /path/to/netmhc-cache/NetMHCpan-4.2-pseudoseqs.dat \
  --netmhciipan-training-pseudo /path/to/netmhc-cache/NetMHCIIpan-4.3-pseudosequence.2023.dat
```

若严格/PBR 资产已经生成，也可运行 `build_netmhc_pseudo_reference.py --official-pseudo-i ... --official-pseudo-ii ...`，随后运行 `build_fallback_reference.py`；fallback builder 会同步刷新 fallback PBR。脚本不会继承旧 `mhcii_model_input.csv` 的裁剪。Human full 记录结合同版 FASTA/status；AF3/PBR 通过固定 IPD alignment 投影 reviewed UniProt 特征，并要求唯一映射回 parent。MHC-I 34 位及 MHC-II alpha15+beta19 alignment mapping 均有版本。只有当 full、PBR 或 deposited partial 可无歧义覆盖全部接触位点时才生成 pseudo；缺位点不会插补。生成输出与冲突报告稳定排序且具有确定性；`metadata.json` 与 `pseudo_validation.json` 记录的是 source input 哈希，并非每个输出文件的哈希。

## 测试

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

单元测试覆盖 strict/fallback 证据隔离、稀有 HLA/H-2、表达后缀与物种过滤、碰撞安全反查、完整模板要求、PBR 双链与边界及 oracle-gated pseudo 核验。严格 AF3/PBR 链均满足 parent-subsequence invariant；无信号肽的 P14437 另按 reviewed UniProt 胞外 topology 核验。确定性另在 release validation 中通过重复完整重建并比较输出哈希检查。

## 重要边界

- `precursor_full` 表示来源确认覆盖完整 coding protein 的记录，通常包含 signal peptide、extracellular domains、TM 和 cytoplasmic tail；真实 null/截短等位基因也可能是完整记录。它不是 AF3 默认输入。
- AF3 的结构输入可直接使用 `af3-soluble-extracellular-v1`；不要把 34-aa pseudosequence 当作可折叠链。
- PBR 虽是连续结构域，但比推荐的完整 AF3 胞外链短，并省略膜近端稳定结构域（MHC-I 还省略 B2M）。
- `source-partial-v1` 是片段，`*-imputed-v1` 是合成构建；两者都不是官方全长等位基因序列。
- pseudo-only scaffold 只约束 34 个位置，并未覆盖全部 peptide/TCR-facing 多态位点；适合初筛，不适合直接解释等位基因特异机制。
- P/G group 只表达特定结构域的一致性，不能当作整条前体蛋白完全一致。
- NetMHC 官方支持成员、可验证 34mer、可用 full/PBR 是三个不同事实；即使名称在支持清单中，没有安全序列来源时仍返回 `not_found`。
- 源数据的许可与使用条款不因本项目代码许可而改变，见 [`NOTICE.md`](NOTICE.md)。

详细审核、TCR 工具对比和 AF3 五链建议见 [`AUDIT_REPORT.md`](AUDIT_REPORT.md)。

## 许可证与数据条款

本包代码采用 [MIT License](LICENSE)。捆绑及派生参考数据保留上游署名和适用条款；MIT 不会重新许可 IPD-IMGT/HLA、UniProt、DTU 或 PUFFIN 数据。本仓库不分发 NetMHCpan/NetMHCIIpan 可执行程序、模型文件或原始训练归档。
