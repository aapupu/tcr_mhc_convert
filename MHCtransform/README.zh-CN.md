# MHCtransform

[English](README.md) | **简体中文**

`MHCtransform` 是 [`aapupu/tcr_mhc_convert`](https://github.com/aapupu/tcr_mhc_convert) 单一仓库中的独立 Python 包。

`MHCtransform` 是一个可审计、离线运行的人/鼠 MHC 转换项目，覆盖：

- MHC-I 重链与物种匹配的 B2M；
- MHC-II α/β 双链；
- IPD/UniProt deposited protein；其中通过官方 Full 状态或明确 UniProt 完整性检查的记录标为 `precursor_full`；
- 适合 AF3 复合物输入的成熟可溶性胞外链 `af3-soluble-extracellular-v1`；
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
| `netmhcpan-34-v1` | MHC-I 成熟重链参考对齐上的 34 个 peptide-contact 位点；不含 B2M | pan-allele 模型 descriptor，不是结构链 |
| `netmhciipan-34-v1` | MHC-II α 链 15 位 + β 链 19 位，按参考对齐拼接；允许 `X` | pan-allele 模型 descriptor，不是结构链 |

位点坐标和来源哈希保存在 [`metadata.json`](src/mhctransform/data/metadata.json)。MHC-II 不能直接用统一的 precursor 下标切片：信号肽长度和对齐插缺不同，因此运行时使用带版本的已发布 pseudosequence artifact，而不是从全长链临时截取。

后备策略逐级累积且默认关闭：`none` → `partial` → `imputed` → `pseudo-scaffold`。34mer 不可逆：2,953 个人 MHC-I 名称仅对应 1,574 个唯一 pseudo，5,620 个人 MHC-II 名称仅对应 1,160 个。因此 imputed 输出只是兼容的结构 scaffold，不是恢复出的生物学全长序列。

## 当前本地参考库

| 内容 | human I | human II | mouse I | mouse II |
|---|---:|---:|---:|---:|
| deposited protein records | 20,066 | 9,419 | 12 | 17 |
| 默认可用的 complete precursor records | 11,283 | 2,525 | 9 | 14 |
| 严格 AF3-ready extracellular-chain records | 11,004 | 2,456 | 9 | 15 |
| 34-aa records | 2,953 | 5,620 | 7 | 8 |

人 HLA 来自固定的 IPD-IMGT/HLA 3.66.0 / commit `b2055ea`，并与同版 `Allele_status.txt` 的官方 `Full/Partial` 状态联合判定；鼠 MHC 和人/鼠 B2M 来自显式 UniProt accession。AF3 边界由 reviewed UniProt topology 投影到同版 IPD 蛋白比对，且逐条验证为前体的连续子序列。34-aa 科学定义来自 NetMHCpan/NetMHCIIpan 文献，本地分发数据采用 PUFFIN 的公开二级镜像。完整来源、accession、内容 SHA-256 和生成日期均在 metadata 中。

deposited records 包含不完整条目；默认查询只返回官方状态为 `Full` 且没有未知残基的人链，或通过 accession/fragment/长度检查的鼠链。`--include-partial` 才会显式纳入片段。当前鼠库是常用、可溯源等位型集合，不声称穷尽全部 H-2 等位基因。

在捆绑的 8,588 个 pseudo 名称范围内，显式后备构建覆盖率为 100%：human I `2,953/2,953`、human II `5,620/5,620`、mouse I `7/7`、mouse II `8/8`。带 `N/L/S/Q/C/A` 表达后缀的人 HLA 默认仍禁止生成 AF3，除非同时明确给出 `--allow-expression-qualified`。

独立 fallback 资产包含 4,161 条 full/fragment 行和 2,097 条 AF3 行；后者由 67 条 deposited exact ectodomain、1,990 条 partial-guided 构建和 40 条 pseudo-only scaffold 组成。

按逐级累积策略统计 pseudo registry 名称覆盖率（本次审核显式允许了带 `N/L/S/Q/C/A` 后缀的 AF3 目标）：

| 范围 | 名称数 | Strict full | `partial` full | `imputed` full | `pseudo-scaffold` full |
|---|---:|---:|---:|---:|---:|
| human I | 2,953 | 1,491 | 2,921 | 2,921 | 2,953 |
| human II | 5,620 | 2,454 | 5,619 | 5,619 | 5,620 |
| mouse I | 7 | 6 | 7 | 7 | 7 |
| mouse II | 8 | 5 | 8 | 8 | 8 |

| 范围 | 名称数 | Strict AF3 | `partial` AF3 | `imputed` AF3 | `pseudo-scaffold` AF3 |
|---|---:|---:|---:|---:|---:|
| human I | 2,953 | 1,491 | 1,551 | 2,921 | 2,953 |
| human II | 5,620 | 2,393 | 2,482 | 5,507 | 5,620 |
| mouse I | 7 | 6 | 6 | 7 | 7 |
| mouse II | 8 | 6 | 6 | 8 | 8 |

AF3 的 `partial` 只接受覆盖完整所选胞外域的 deposited 序列，绝不会给短片段自动补尾；`imputed` 要求有可用 deposited 引导并排除 pseudo-only 构建。只有显式选择最低证据级别 `pseudo-scaffold` 才覆盖完整 34-aa registry。

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
    "HLA-B*08:03",
    mhc_class="I",
    definition_id=AF3_DEFINITION,
    fallback="imputed",
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

# 仅从已有严格资产重建确定性的 fallback 数据
python3 tools/build_fallback_reference.py
```

脚本不会从旧版 `mhcii_model_input.csv` 继承裁剪结果。人 HLA full registry 用同版 FASTA 和 `Allele_status.txt` 联合构建；AF3 链用固定版 IPD protein alignment 把 exact allele 映射到 reviewed UniProt extracellular topology，并要求结果可唯一定位回 parent precursor；pseudo artifact 再用官方 2009→现行命名表和 deleted/renamed 表规范旧名。鼠链按显式 accession 获取；可变的 UniProt REST 批次另外用预期 SHA-256 固定，内容变化会 fail closed，不能静默改变边界。所有生成文件均稳定排序并写入来源哈希。二级 pseudo 镜像中的 5 个源冲突没有静默保留或任选首条，而是按可审计规则写入 `quarantined_pseudosequences.csv`。

## 测试

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

41 项测试在原严格行为之外，还覆盖 partial/imputed 证据隔离、exact-ectodomain 提升、稀有 HLA 与 H-2 后备、P14437 精确 I-A-s topology、表达后缀与物种过滤、不会漏掉 strict 序列碰撞的显式后备反查，以及全部 8,588 个 pseudo 名称同时具有 AF3/full 后备、所有 imputed full 都使用完整 precursor 模板这两项 invariant。严格 AF3 链继续满足 parent-subsequence invariant；无信号肽的 P14437 另外按其精确 UniProt 胞外 topology 验证。

## 重要边界

- `precursor_full` 表示来源确认覆盖完整 coding protein 的记录，通常包含 signal peptide、extracellular domains、TM 和 cytoplasmic tail；真实 null/截短等位基因也可能是完整记录。它不是 AF3 默认输入。
- AF3 的结构输入可直接使用 `af3-soluble-extracellular-v1`；不要把 34-aa pseudosequence 当作可折叠链。
- `source-partial-v1` 是片段，`*-imputed-v1` 是合成构建；两者都不是官方全长等位基因序列。
- pseudo-only scaffold 只约束 34 个位置，并未覆盖全部 peptide/TCR-facing 多态位点；适合初筛，不适合直接解释等位基因特异机制。
- P/G group 只表达特定结构域的一致性，不能当作整条前体蛋白完全一致。
- 当前 pseudosequence artifact 与最新 IPD 全长快照覆盖范围不同；没有 exact pseudo 记录时返回 `not_found`。
- 源数据的许可与使用条款不因本项目代码许可而改变，见 [`NOTICE.md`](NOTICE.md)。

详细审核、TCR 工具对比和 AF3 五链建议见 [`AUDIT_REPORT.md`](AUDIT_REPORT.md)。

## 许可证与数据条款

本包代码采用 [MIT License](LICENSE)。捆绑及派生的参考数据仍保留上游署名和适用条款；MIT 许可证不会重新许可 IPD-IMGT/HLA、UniProt 或 PUFFIN 派生数据。
