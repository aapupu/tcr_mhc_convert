# MHC 转换审核、TCR 工具对比与 AF3 五链输入建议

首次审核：2026-11-03；fallback 扩展复核：2026-11-20；DTU pseudo 主来源与 PBR 硬覆盖复核：2026-11-21

## 结论

旧实现不能作为“人/鼠 MHC-I/MHC-II、全长与 34-aa、双向转换”的正确实现继续使用。它的人 HLA-I 34-aa 数值抽查是对的，但完整范围、命名、反向语义和来源版本不足；MHC-I 默认链不是全长；MHC-II 裁剪表存在系统性错误；mouse 与 MHC-II pseudo 基本缺失。

本轮复核后的新 `MHCtransform` 在其明确声明的覆盖范围内可以通过审核：固定来源、保留歧义、strict full/PBR/pseudo/AF3 与 fallback 分层，支持人/鼠 MHC-I/MHC-II 正反向查询。它不是全 IPD-MHC 物种全集；mouse 是显式 accession 的常用集合。Pseudo 名称空间已对齐 NetMHCpan 4.2f 与 NetMHCIIpan 4.3k 的公开支持清单，但“官方支持名称”不等于“已有可安全分发的序列”。医学/移植决策仍应回到当前官方数据库验证。

新实现将以下概念彻底分开：

1. `precursor_full`：IPD 官方 sequence status 为 `Full` 且无未知残基，或经显式 UniProt accession/fragment 检查确认覆盖完整 coding protein 的记录；其余 deposited protein 仅在 opt-in 时返回；
2. `af3-soluble-extracellular-v1`：成熟、完整的可溶性胞外链，去 signal/TM/cytoplasmic，但保留膜近端 connector；
3. 连续 PBR：MHC-I 为一条 alpha1+alpha2；MHC-II 为独立 alpha1 与 beta1 两条链。它用于肽结合槽比较与 pseudo 核验，不替代完整胞外 AF3 链；
4. 有 `definition_id` 的 34-aa descriptor：用于泛等位模型，不是结构链；
5. 显式 fallback：`source-partial-v1` 是 deposited fragment；`precursor-imputed-v1` 与 `af3-soluble-extracellular-imputed-v1` 是保留 deposited/pseudo 约束、其余位置由同物种同位点完整模板提供的合成 scaffold。它们不进入 strict 索引，也不是官方等位基因全长序列。

正向查询不任意挑选 first row。严格反查返回全部 `exact_*` candidates；fallback 构建默认不参与反查，显式纳入后只返回 `compatible_*`，防止把合成序列误写成生物学精确匹配。

### 本轮发现并修复的高风险问题

- 旧 canonicalization 曾把 `A*03:03N` 的 `N` 后缀剥掉，导致不存在的裸 `HLA-A*03:03` 命中 null allele；另一种不安全简化会把官方高字段 `A*02:01:01:02L` 人工变成并不存在的 `A*02:01L`。现在 `N/L/S/Q/C/A` 连同完整官方字段全部保留，AF3 默认排除这些表达状态；裸名与带后缀名的索引严格隔离。
- full lookup 曾把任意不存在的高字段名自动投到两字段，例如 `A*02:01:999:999`。现在默认 `not_found`；只有显式 `allow_two_field_projection` 才投影，且结果状态为 `projected_*`。如果精确名称是官方已知的 Partial（例如 `A*02:01:17`），AF3 即使开启投影也 fail-closed，不能冒充同两字段的完整链。Pseudo 的高字段投影也先验证该 deposited name 确实存在。
- AF3 builder 曾有两个潜在错误：IPD 可视化 alignment 的 partial/stop-loss 行不能按普通空白切割；UniProt topology 也不能只按长度套到 alignment 第一行。现在 parser 保留对齐列，且每个 locus 必须先找到与 UniProt topology reference **全序列完全相同**的 IPD row，再向 exact deposited allele 投影。
- AF3 provenance 现在一行对应一个 parent full record 和一个明确 `start-end`，避免 parent ID 列表与 span 列表失配。

本轮扩展后保留 29,514 个 deposited protein records、13,484 个 strict AF3-ready chains、13,482 个 strict PBR chains 与 24,465 个 pseudo records（MHC-I 15,489；MHC-II 8,976）。严格 AF3/PBR 都验证为 parent 的连续子序列；无信号肽的 mouse I-A-s alpha 链 P14437 另按 reviewed UniProt 胞外 topology 验证。Pseudo 输出另有逐条 evidence、official-support、validation status、来源记录与质量标志。

## PBR 与 34-aa pseudo 的关系

新增 PBR 是“完整肽结合槽结构域”，而不是另一套 pseudo：

- MHC-I 输出一条连续 heavy-chain alpha1+alpha2，经典 HLA/H-2 通常 182 aa，不含 alpha3 与 B2M；
- MHC-II 输出两条独立链：alpha1（人类常见 84–96 aa）与 beta1（90–95 aa），绝不人工拼接；
- NetMHCpan-I 的 34 个接触位点全部位于 alpha1+alpha2；NetMHCIIpan-II 的 alpha15/beta19 全部位于各自 alpha1/beta1，因此 PBR 是提取与反向核验 pseudo 的合适连续上位表示；
- 构建器要求每个 pseudo 组件的全部接触位点都能投影到一个 PBR；所有非 `X` 残基必须逐位相同，DTU `X` 只作为 alignment/unknown 通配符记录，绝不作为氨基酸写入蛋白链；
- pMHC–TCR 的 AF3 默认仍使用完整可溶性胞外链。仅用 PBR 会省略膜近端稳定域，并可能暴露非天然截面；MHC-I 还会缺失 B2M。

人类边界采用 reviewed UniProt domain REGION，并通过同位点 IPD alignment 投影到具体等位基因：MHC-I 182；DRA 84；DQA1 96（少数 indel 为 95）；DPA1 87；DRB1/3/4/5 95；DQB1 94；DPB1 92 aa。小鼠按具体 accession/位点投影，不把一个固定长度套到所有 H-2。

## NetMHCpan 4.2e pseudo、4.2f support 与 NetMHCIIpan 4.3k 对齐

项目提交了支持名称快照、用户提供的 DTU pseudosequence 文本及 SHA-256；不提交 DTU 可执行程序、模型文件或完整训练归档。NetMHCpan 的 `4.2e` 是 supplied pseudo artifact 的版本，`4.2f` 是独立 support list 的版本，两者不能混写：

- NetMHCpan 4.2f support：人/鼠审核范围 10,395 个支持名称，本地覆盖其中 8,744 个；
- NetMHCIIpan 4.3k support：DR、DPA1×DPB1、DQA1×DQB1 及小鼠共 11,344 个支持分子，本地覆盖其中 8,976 个；
- 扩展 registry 共 24,465 行：DTU 直接 pseudo 5,848、DTU component composition 3,339、`derived_exact_pbr` 10,808、`derived_contact_complete_partial` 3,844、`legacy_secondary` 626；组件 pseudo 4,390 条，另有 344 条冲突/不可接受记录隔离。

DTU training pseudo 现在是直接主来源：216 条人/鼠 MHC-I 与 5,632 条人/鼠 MHC-II 行以 `dtu_official_pseudosequence` 证据进入 registry；由 DTU 认证的 MHC-II 链组件还能组合出 3,339 行。Deposited PBR/full/partial 派生与 legacy 行继续提供更广名称空间，但不能覆盖 DTU 直接值。`H-2-Kq` 因此不再缺失。DQ/DP 的笛卡尔支持空间只代表模型允许的链组合，不代表 haplotype 或体内稳定性证据。

命名和冲突不静默处理：历史 `HLA-DQA1*01:07` 与 `HLA-DRB1*15:13` 分别解析为当前 `HLA-DQA1*01:07Q` 与 `HLA-DRB1*15:13Q`；compact `DRB1_13100` 解析为 `HLA-DRB1*13:100` 而非 `*131:00`。DTU pseudo 对 pseudo 表示具有优先级，但 current-IPD strict 蛋白/PBR 不会因此被篡改；冲突时生成带 synthetic/provenance 标记的兼容 PBR fallback，并且只有显式 `--require-pseudo-match` 才选用。

`pseudo_components.csv` 的 4,390 行是链组分审计表。4 个争议 DPB1 beta 组件按 DTU 直接值进入 pseudo registry；它们与 current deposited derivation 的冲突保留在质量标志与 quarantine 中。Pair 级 pseudo 与 provenance 均保留，但仍不能据此宣称某个 DP/DQ 组合具有生物学 haplotype 证据。

`pseudo_pbr_projection.csv` 为每个唯一 pseudo chain component 保存选定 PBR、全部接触坐标、投影残基、`X` wildcard 位置及 strict/fallback 来源。构建时未覆盖组件数必须为 0；运行时 `allele-to-pbr`/`molecule-to-pbr --require-pseudo-match` 只返回该审计表认证的 PBR。MHC-II 的 alpha1 与 beta1 分别认证，任一链未通过都不能把整个分子报告为 pseudo-compatible。

## 基于 34-aa pseudo 的全面 fallback 复核

扩展后的 34-aa registry 仍不能把 pseudo 唯一恢复成 full/PBR/AF3；多个名称可共享同一 34mer，且 pseudo 不定义接触位点以外的序列。

Fallback 能提供可追踪的结构 scaffold，但不表示从 34 aa 恢复了真实全长。本次固定构建有 full/AF3/PBR fallback 13,219/6,628/6,628 条；18,080 个唯一 pseudo 组件全部具有 PBR 证书，其中 strict 11,559、fallback 6,521、未覆盖 0。坐标审计仅允许唯一且最大 1 aa 的恢复：17 个组件发生恢复，其中 12 个含 alignment-gap contact。Exact ectodomain 表示 deposited fragment 覆盖所选 AF3 全区，也不代表它同时覆盖 signal/TM/cytoplasmic precursor；机器可读的当前计数仍以 `metadata.json` 为准。

多个等位基因可共享同一 34mer，且这些位置不能确定非 pseudo 区域或全部 TCR 接触面。因此 fallback 采用累积、默认关闭的四级策略：

1. `none`：只返回 strict 官方/可核验序列；
2. `partial`：仅返回 deposited fragment，绝不标为 full；
3. `imputed`：有 fragment 时优先保留其对齐残基，再用同物种、同位点、同链型完整模板补齐；
4. `pseudo-scaffold`：没有 deposited sequence 时，仅用非 `X` pseudo 位点约束同位点模板，是最低证据级别。

所有构建都记录 `evidence_level`、来源记录、模板、observed/imputed 区间、pseudo 距离、未确定的 `X` 位点及质量标志。带 `N/L/S/Q/C/A` 表达后缀的人 HLA 默认仍不能进入 AF3 fallback，除非使用者额外显式确认。

截图中的例子经逐项复核：`HLA-B*08:03`、`HLA-B*45:06`、`HLA-A*02:50`、`HLA-A*03:19` 与 `H-2-Lq` 均使用 deposited partial 引导的 AF3 构建；`H-2-IAs` 的 α 链 P14437 1–195 是 exact extracellular sequence，不需要补全；`H-2-IAq` 与 `H-2-IAu` 仅 α 链需要 imputed，β 链为 exact。合成 full 的模板 parent 必须通过 strict precursor completeness 检查，不能把 P14437 这类 signal-free mature membrane record 冒充 full precursor。

## 旧版 MHC 部分审核

审核对象是保留在 legacy 目录的 [`convert_tcr_pmhc_full_sequences.py`](../legacy/tcr_design_v1/convert_tcr_pmhc_full_sequences.py) 及其默认数据表。

### 1. `--species mouse` 没有切换 MHC 数据

`--species` 只选择 TCR V/J/C reference。默认 MHC-I 表有 22,780 行、MHC-II 表有 8,757 行，但两表中 `H-2` 记录均为 0。因此旧命令即使指定 mouse，也不能完成 mouse MHC 转换。

### 2. 旧版所谓 `full_aa` 不是同一个生物学概念

- HLA-A*02:01:01:01 的默认 MHC-I `full_aa` 是 274 aa 的成熟可溶性 ectodomain；同一行 `full_aa_original` 才是 365 aa precursor。
- B2M 默认是 99 aa mature chain，而不是 119 aa precursor。
- 变量名和 CLI 帮助仍称其为 `full_aa`，容易把结构建模链误当成表达意义全长。

### 3. MHC-II “AF-ready” 表不可靠

在 `mhcii_model_input.csv` 的 8,757 行中：

- 270 行 `full_aa` 少于 150 aa；最短仅 99 aa；
- 3,214 行以 `M` 开头；起始 `M` 本身不能证明存在 signal peptide，但代表例 `HLA-DRA*01:01:01:01` 确实保留明显疏水 signal peptide `MAISGVPVLGFFIIAVLMS...`；
- 多个 DRB 条目仅剩从结构域中部开始的 105 aa 片段，却仍被标记 `model_ready=True`。

所以不能通过改名或小修补继续使用这张裁剪表。

### 4. two-field fallback 会隐藏歧义

旧 `AlleleLookup.add()` 对 two-field key 使用 `setdefault`，只保留读到的首条记录。两字段名称可能对应多个 deposited precursor；34-aa descriptor 更经常由许多 allele 共享。任意选第一条会把数据顺序伪装成生物学唯一性。

### 5. pseudo 覆盖不完整且没有反向

- `doc/pseudo_HLAI.csv` 是 human HLA-I-only；A*02:01 等抽查值与 NetMHCpan 34 位定义一致。
- 没有 mouse MHC-I、human/mouse MHC-II 34-aa 实现。
- 没有 full sequence → allele 或 pseudo → allele。
- 标准 `HLA-A*02:01` 与旧表的无星号命名并不总能直接匹配。
- `script/dataprocess.py` 的 HLA-II pseudo 路径指向一个不存在的文件。

## 新实现的科学定义

### MHC-I 34 aa

采用 NetMHCpan 类型的成熟重链参考位置（1-based）：

`7, 9, 24, 45, 59, 62, 63, 66, 67, 69, 70, 73, 74, 76, 77, 80, 81, 84, 95, 97, 99, 114, 116, 118, 143, 147, 150, 152, 156, 158, 159, 163, 167, 171`

定义只描述 heavy chain，不包含 B2M。这些也是成熟重链的参考对齐位置，不是 precursor 的通用数组下标；例如人 HLA 与鼠 H-2 的 signal peptide 长度不同。原始 34 位由经典 HLA-A/B 结构与 A/B/C 多态位点建立；用于 HLA-E/F/G 或 mouse Qa/T 等非经典分子应视为外推，而不是同等验证。原始方法见 [Nielsen et al., PLoS ONE 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC1949492/)。

### MHC-II 34 aa

采用 NetMHCIIpan 类型的 α15 + β19：

- α：`9, 11, 22, 24, 31, 52, 53, 58, 59, 61, 65, 66, 68, 72, 73`
- β：`9, 11, 13, 26, 28, 30, 47, 57, 67, 70, 71, 74, 77, 78, 81, 85, 86, 89, 90`

这些是成熟链的**参考对齐列**，不是每条 precursor 可直接使用的 Python index。NetMHCIIpan artifact 用 `X` 表示参考对齐中的 gap/unknown token；结构性 deletion 与数据库序列缺失不能混为一谈。DRB shorthand 历史上隐含 DRA*01:01；DQ/DP 必须显式给 α+β。artifact 中出现一个 DQA–DQB 或 DPA–DPB pair 只表示模型序列对，不证明该组合具有 haplotype 或体内稳定性证据。方法依据见 [Karosiene et al., Immunogenetics 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3809066/)。旧版 NetMHCIIpan 的 β-only 21-aa 定义不能与此 34-aa 定义混用。

### 数据与命名

- human deposited protein 与官方 Full/Partial 状态：固定版 [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) `hla_prot.fasta` + `Allele_status.txt`；IPD 的 [FAQ](https://www.ebi.ac.uk/ipd/imgt/hla/about/help/faq/) 明确说明很多命名等位基因只有 partial CDS；
- human 等位基因命名规则： [WHO HLA nomenclature](https://hla.alleles.org/pages/nomenclature/naming_alleles/)；
- mouse gene/allele/strain 命名以 MGI 为权威层级，MHC 名称/蛋白映射参考 [IMGT mouse MHC nomenclature](https://www.imgt.org/IMGTrepertoireMH/LocusGenes/nomenclatures/mouse/MHC/Mu_MHCnom.html)；
- mouse/B2M 蛋白：显式 [UniProt](https://www.uniprot.org/) accession；
- 34-aa 标准来源是上述 NetMHCpan/NetMHCIIpan 论文；支持名称使用 DTU 4.2f/4.3k 公开快照；随库保存的 NetMHCpan 4.2e training pseudo 与 NetMHCIIpan 4.3k `pseudosequence.2023.dat` 是直接主来源。PBR/full/partial 精确派生及 [PUFFIN public repository](https://github.com/gifford-lab/PUFFIN) legacy secondary 用于扩展名称空间。生成输出稳定排序且具有确定性；metadata/validation report 记录 source input 哈希，而不是宣称每个输出文件都有记录哈希。

MHCflurry 是 MHC-I-only，并同时存在 34/37/39 位 artifact；34 位是 NetMHCpan-derived 表示，37/39 位是其他模型表示。模型权重依赖精确定义，不能只看“pseudo”字段名就互换。其实现说明见 [MHCflurry pseudosequences module](https://openvax.github.io/mhcflurry/_modules/mhcflurry/pseudosequences.html)。

Pseudo 构建器保留 HLA 冒号字段边界，避免 `A*02:11:02` 与 `A*02:1102` 碰撞；所有来源冲突及不完整接触位点都写入 quarantine/validation report，不静默任选首条。MHC-II molecule resolver 只允许 DRalpha–DRbeta、DQalpha–DQbeta、DPalpha–DPbeta 及对应 mouse I-A/I-E 家族配对。

## AF3 链转换的审核结果

AF3 官方输入格式只要求每个 protein entity 有独立 ID 和单字母氨基酸序列，并没有发布 pMHC–TCR 的统一裁剪坐标。因此 `af3-soluble-extracellular-v1` 是本项目的工程约定，不应写成“AF3 官方序列”。其默认策略是保留完整成熟胞外折叠单元，避免裁掉原有结构域后暴露新的假 docking surface。

| 例子 | precursor | AF3 边界（1-based） | 长度 | 下一段证据 |
|---|---:|---:|---:|---|
| HLA-A*02:01 heavy | 365 | 25–308 | 284 | 309–332 为 HLA-A TM |
| human B2M | 119 | 21–119 | 99 | signal 1–20 |
| H-2-Kb heavy | 369 | 22–305 | 284 | 306 起为疏水 TM |
| mouse B2M | 119 | 21–119 | 99 | signal 1–20 |
| HLA-DRA*01:01 | 254 | 26–216 | 191 | 217 起为 TM |
| HLA-DRB1*15:01 | 266 | 30–227 | 198 | 228–248 为 TM |
| H-2-IAb α | 256 | 24–218 | 195 | 219 起为 TM |
| H-2-IAb β | 265 | 28–226 | 199 | 227 起为 TM |

HLA-A 的 topology 来源是 reviewed UniProt locus record P04439，但 exact A*02:01 sequence 来自 IPD `HLA00005`；当前 P04439 canonical sequence 不是 A*02:01，不能把该 UniProt sequence 直接冒充目标 allele。实现用 exact sequence-matched IPD alignment row 投影 25–308，并再次确认输出是目标 precursor 的唯一连续子串。

示例：

- HLA-A*02:01 AF3：`GSHSMRYFFTSVSRPGRGEP...GLPKPLTLRWEPSSQPTIPI`（284 aa）；MHC-I molecule 同时输出 human B2M `IQRTPKIQVYSRHPAENGKS...CRVNHVTLSQPKIVKWDRDM`（99 aa）。
- HLA-DRA*01:01 / DRB1*15:01：分别输出 191/198 aa 两条链，不能拼接。
- H-2-IAb：分别输出 H2-Aa*b 195 aa 与 H2-Ab1*b 199 aa 两条链。

旧 274-aa HLA-I chain 截至约 precursor 298，少了 299–308 的 membrane-proximal connector；它比 precursor 更接近结构输入，但不是本项目的新默认。AF3 reverse lookup 仍是 many-to-many：多个 deposited allele 可共享完全相同的 extracellular chain，程序返回全部 exact candidates。

## TCRconvert、旧 converter、TCRtransform、Stitchr

首先要消歧名称：公开的 [seshadrilab/TCRconvert](https://github.com/seshadrilab/tcrconvert) 只转换 10x、Adaptive、IMGT 等 V/D/J/C **基因命名**，不重建蛋白，也不存在“保留哪些 CDR/FR”的问题。真正可与本地 `TCRtransform` 比较的公开重建器是 [Stitchr](https://github.com/JamieHeather/stitchr)。

| 能力/区域 | public TCRconvert | 本地旧 converter | `TCRtransform` | Stitchr |
|---|---|---|---|---|
| 主要任务 | gene name 映射 | TCR AA 重建 + MHC 查表，输出 A–E 五链 | 结构用 TCR AA 重建与反查 | 等位基因感知的 nucleotide/cDNA stitching |
| leader/signal peptide | 不适用 | 不保留 | 不保留 | 默认保留；可 `--no_leader` |
| V 区 | 不适用 | 成熟 V 起点到重排边界 | 保留 FR1–FR3、CDR1、CDR2、CDR2.5；按 junction 裁 V 末端 | 等位基因级 V，在 nucleotide 层裁剪 |
| CDR3 | 不适用 | inclusive-anchor junction | 从 V 侧 C104 到 J 侧 F/W118 的 inclusive-anchor `junction_aa`，只保留一次 | 接受同类 junction 并完成重排 |
| D 区 | 不适用 | 位于输入 CDR3 内 | β 链 D/N/P 信息保留在输入 CDR3，不单独猜 TRBD | 在 junction 中处理 |
| J/FR4 | 不适用 | 保留 | 保留锚点后的 J tail / FR4 | 保留并完成 J–C splice |
| constant | 不适用 | 结构用 C-domain | α/β 的 extracellular C-domain；α 补 TRAJ–TRAC split-codon residue；β 按 TRBJ1/2 选 TRBC1/2 | 完整 constant coding exons |
| connecting/TM/cytoplasmic | 不适用 | 不保留 | 不保留 | 默认完整 cDNA 中保留 |
| 输出 | gene labels | AA 五链 | AA TCR 链与区域注释 | nucleotide，亦可翻译 AA |

因此 `TCRtransform` 文档中的“full chain”应理解为**成熟、结构建模用 V+C ectodomain**，不是含 signal peptide、TM 和胞内尾的表达载体全长。

CDR1/CDR2 使用 IMGT 定义；CDR2.5 是 tcrdist 的 pMHC-facing/HV4-like loop（IMGT 对齐列 81–86），不是第四个经典 IMGT CDR。严格 CDR3-IMGT 不含 C104 与 F/W118 锚点，而本项目字段实际是常见 repertoire 意义的 inclusive-anchor junction。

现有 [`TCRtransform/validation/REPORT.md`](../TCRtransform/validation/REPORT.md) 的 400 链测试显示：可比的 CDR1/2/2.5 为 399/399 一致；双方成功时 CDR3 359/359 只出现一次；共享模型区 297/359（82.7%）逐残基一致。主要差异来自 `TCRtransform` 的 gene-level reference fallback 丢失等位基因特异序列，以及默认重叠阈值更严格。

## AlphaFold 3 的 pMHC–TCR 五链输入区域

### 先区分“官方要求”和“工程建议”

[AF3 input specification](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md) 要求把不同蛋白作为独立 entity，但没有给 TCR/pMHC 专用裁剪坐标。因此下述是基于免疫复合物结构与建模研究的默认建议，不是 Google 官方规定。

### 推荐的五条链

| 链 | MHC-I 复合物 | MHC-II 复合物 | 推荐区域 |
|---|---|---|---|
| 1 | TCR α | TCR α | 成熟 extracellular Vα + Cα，按具体 accession 的 extracellular/topological annotation 保留至 TM 前；包含 FR1–FR4、全部 CDR 与 constant Ig domain；去 signal/TM/cytoplasmic |
| 2 | TCR β | TCR β | 成熟 extracellular Vβ + Cβ；同上 |
| 3 | MHC-I heavy | MHC-II α | I：按具体 accession 保留完整 extracellular α1+α2+α3 至 TM 前；II：完整 extracellular α1+α2 至 TM 前；去 signal/TM/cytoplasmic |
| 4 | mature B2M | MHC-II β | I：完整 mature B2M；II：按具体 accession 保留完整 extracellular β1+β2 至 TM 前；去 signal/TM/cytoplasmic |
| 5 | peptide | peptide | 使用实验/生物学上的完整呈递肽；MHC-II 已知天然 flanks 时保留，不要只输入 9-aa core |

表中链号只是项目约定，AF3 不要求这一顺序；截取边界必须逐 accession 确认，没有通用的“短 stalk”坐标。

不要把非连续的 34-aa descriptor 作为 AF3 蛋白链；它没有可折叠的主链拓扑。

默认也不建议只给 TCR V domains 或只给 MHC peptide-binding domains。[TCRdock/eLife 2023](https://elifesciences.org/articles/82813) 作者指出，只保留直接接触域虽然更快，但可能在缺失域原有界面处制造 decoy docking surface；因此这里把完整可溶性 ectodomain 作为保守默认。该证据来自 AlphaFold-Multimer，应视为对 AF3 的工程外推，而不是 AF3 已验证的硬规则。近期 AF3 工作直接支持把 peptide、两条 MHC 链和 TCRα/β 作为五个 constituent chains 建模，见 [Science Advances 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13252052/)，但不规定统一截取坐标。

### 与本项目输出的关系

- `TCRtransform` 的默认输出接近推荐的 TCR V+C 结构链，可作为起点；它可能不含 connecting/stalk，仍需按具体 accession 检查链尾边界。
- `MHCtransform precursor_full` **不能原样作为默认 AF3 链**，因为它包含 signal peptide、TM 和 cytoplasmic tail；应改用新增的 `af3-soluble-extracellular-v1`。
- 旧 MHC-I 的 274-aa soluble chain 比 precursor 更接近 AF3 方向，但少了 HLA-A 299–308 的近膜 connector；旧 MHC-II 自动裁剪同样不可信。新 AF3 registry 已按 exact IPD alignment + UniProt topology 重建，不继承这些旧裁剪。
- `molecule-to-af3` 对 MHC-I 输出 heavy + 同物种 mature B2M 99 aa；对 MHC-II 输出 α+β 两条链。

### AF3 结果检查

至少运行多个 seed，并重点看。Science Advances 2026 的具体配置是 MSA/templates 开启、5 个 seeds、每 seed 1 个 diffusion sample；这是该研究配置，不是 Google 官方固定最佳实践。

- peptide–MHC 与 TCR–pMHC 的 pairwise ipTM / inter-chain PAE；
- peptide 的 pLDDT 与相对 MHC groove 的 PAE；
- 不同 seed 是否给出一致的 docking orientation；
- 是否出现链尾、裁剪面或假性暴露疏水面主导的界面。

只看全局 pLDDT 不足以判断五链界面是否可信。AF3 置信度解释可参考 [EMBL-EBI AlphaFold 3 quality guide](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/how-to-assess-the-quality-of-alphafold-3-predictions/)。

