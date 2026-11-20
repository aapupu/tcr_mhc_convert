# tcr_mhc_convert

[English](README.md) | **简体中文**

`tcr_mhc_convert` 是一个可审计、可离线运行的 Python 单一仓库，用于在人和小鼠的 MHC/TCR 注释与蛋白质序列之间进行转换；反向转换仅在能够精确匹配时返回结果。

本仓库由原 `aapupu/tcr_design` 重命名而来。

## 两个独立工具

| 工具 | 主要功能 | 命令行 |
|---|---|---|
| [`MHCtransform`](MHCtransform/) | 人/鼠 MHC-I、MHC-II 名称与前体全长、AF3 胞外链、连续肽结合区（PBR）及 34-aa pseudosequence 的严格双向转换，并提供与严格库隔离、需显式开启的后备；MHC-I 自动配物种匹配 B2M，MHC-II 输出 alpha/beta 双链 | `mhctransform` |
| [`TCRtransform`](TCRtransform/) | 人/鼠 alpha/beta TCR 的 V + CDR3 + J + C 胞外链重建、CDR1/2/2.5/3 注释、CDR3 边界解析和完整链反查 | `tcrtransform` |

两个包运行时都只依赖 Python 标准库，支持 Python 3.9 及以上版本。它们共享一个 GitHub 仓库，但保留独立的 Python 包名、版本、API 和命令行，因此可以单独或一起安装。

## 安装

```bash
git clone https://github.com/aapupu/tcr_mhc_convert.git
cd tcr_mhc_convert

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ./MHCtransform -e ./TCRtransform
```

## 快速示例

### MHC allele -> 前体全长、肽结合区、AF3 链或 34-aa pseudo

```bash
# HLA-A 重链前体
mhctransform allele-to-full 'HLA-A*02:01' --class I

# 适合 AF3 的 HLA-A 胞外重链 + human B2M
mhctransform molecule-to-af3 'HLA-A*02:01' --class I

# 适合 AF3 的 MHC-II alpha/beta 双链
mhctransform molecule-to-af3 \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# 连续肽结合结构域：MHC-I 输出一条 alpha1+alpha2 链
mhctransform molecule-to-pbr 'HLA-A*02:01' --class I

# MHC-II PBR 输出独立的 alpha1 与 beta1 两条链
mhctransform molecule-to-pbr \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# 要求 PBR 的接触位点逐位复现 DTU pseudo
mhctransform molecule-to-pbr 'HLA-A*02:01' \
  --class I --require-pseudo-match

# 稀有等位基因：保留 deposited 残基，只补齐 AF3 缺失位置
mhctransform molecule-to-af3 'HLA-B*08:03' \
  --class I --fallback imputed

# NetMHCpan 兼容的 34-aa descriptor
mhctransform allele-to-pseudo 'HLA-A*02:01' \
  --class I --definition netmhcpan-34-v1
```

反向查询使用 `full-to-allele`、`af3-to-allele`、`pbr-to-allele` 和 `pseudo-to-allele`。严格序列会返回全部 exact 候选；后备构建默认不参与反查，只有显式开启时才返回 `compatible_*` 候选，且不声称生物学精确匹配。详细定义、覆盖范围和审核结论见 [`MHCtransform/README.zh-CN.md`](MHCtransform/README.zh-CN.md) 与 [`MHCtransform/AUDIT_REPORT.md`](MHCtransform/AUDIT_REPORT.md)。

### TCR V/J/CDR3 -> 适合 AF3 的 alpha/beta 胞外链

准备至少含有 `TRAV, TRAJ, CDR3A, TRBV, TRBJ, CDR3B` 的 CSV：

```bash
tcrtransform build TCRtransform/examples/tcr_input.csv \
  --species-column species \
  -o tcr_built.csv
```

只注释 CDR 区域可用 `tcrtransform annotate`；从已有完整胞外链反查 V/J/C 和 CDR 可用 `tcrtransform parse`。详细说明见 [`TCRtransform/README.zh-CN.md`](TCRtransform/README.zh-CN.md)，对照 Stitchr/tcrdist3 的验证见 [`TCRtransform/validation/REPORT.md`](TCRtransform/validation/REPORT.md)。

## AlphaFold 3 的 pMHC-TCR 五链输入

建议把每条分子作为独立 protein entity 输入，不拼接 linker，也不要把 34-aa pseudo 当作结构链：

| AF3 entity | MHC-I | MHC-II |
|---|---|---|
| 1 | peptide | peptide |
| 2 | MHC 重链胞外域 | MHC alpha 链胞外域 |
| 3 | B2M | MHC beta 链胞外域 |
| 4 | TCR alpha 成熟 V+C 胞外链 | TCR alpha 成熟 V+C 胞外链 |
| 5 | TCR beta 成熟 V+C 胞外链 | TCR beta 成熟 V+C 胞外链 |

`MHCtransform` 的 `af3-soluble-extracellular-v1` 会移除 signal peptide、跨膜区和胞内尾，保留全部胞外折叠域及膜近端连接区；`TCRtransform` 输出成熟的 TCR V+C 胞外结构链。AF3 官方输入格式只要求分别提供各条氨基酸链，并未规定统一裁剪坐标；上述范围是本仓库经过来源核对和测试的工程 profile。Imputed 与 pseudo-scaffold MHC 链可用于明确标注前提下的结构初筛，但不能解释为未知等位基因的真实全长序列。

## 测试

```bash
PYTHONPATH=MHCtransform/src:TCRtransform/src \
  python -m unittest discover -s MHCtransform/tests -p 'test*.py' -v

PYTHONPATH=MHCtransform/src:TCRtransform/src \
  python -m unittest discover -s TCRtransform/tests -p 'test*.py' -v
```

根目录的 [GitHub Actions workflow](.github/workflows/tests.yml) 会在 Python 3.9 和 3.11 上安装并测试两个包。

## 仓库结构

```text
tcr_mhc_convert/
|-- MHCtransform/          # mhctransform 包、数据、审核和测试
|-- TCRtransform/          # tcrtransform 包、验证和测试
|-- legacy/tcr_design_v1/ # 为追溯保留的旧版文件
|-- .github/workflows/     # 两个包的持续集成
|-- LICENSE                # MIT 代码许可证
|-- README.zh-CN.md        # 本页
`-- README.md              # 默认英文首页
```

## 许可证与数据条款

本仓库代码采用 [MIT License](LICENSE)。每个子项目的 `NOTICE.md` 和数据目录元数据记录了 IPD-IMGT/HLA、UniProt、tcrdist3、NetMHCpan 4.2f/NetMHCIIpan 4.3k 公开支持清单、随库保存的 NetMHCpan 4.2e/NetMHCIIpan 4.3k pseudosequence 快照及旧 PUFFIN 镜像等来源、固定版本和校验值。DTU pseudosequence 文本是 MHC registry 的直接主输入；仓库不分发 DTU 可执行程序、模型文件或完整训练归档。MIT 许可证只适用于本仓库代码，不会重新许可所捆绑或派生的生物参考数据。再分发数据前请核对各上游条款。
