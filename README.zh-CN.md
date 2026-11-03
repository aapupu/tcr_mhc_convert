# tcr_mhc_convert

[English](README.md) | **简体中文**

`tcr_mhc_convert` 是一个可审计、可离线运行的 Python 单一仓库，用于在人和小鼠的 MHC/TCR 注释与蛋白质序列之间进行转换；反向转换仅在能够精确匹配时返回结果。

本仓库由原 `aapupu/tcr_design` 重命名而来。

## 两个独立工具

| 工具 | 主要功能 | 命令行 |
|---|---|---|
| [`MHCtransform`](MHCtransform/) | 人/鼠 MHC-I、MHC-II 的 allele/name 与前体全长、AF3 胞外链、34-aa pseudosequence 双向精确转换；MHC-I 自动配物种匹配 B2M，MHC-II 输出 alpha/beta 双链 | `mhctransform` |
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

### MHC allele -> 前体全长、AF3 链或 34-aa pseudo

```bash
# HLA-A 重链前体
mhctransform allele-to-full 'HLA-A*02:01' --class I

# 适合 AF3 的 HLA-A 胞外重链 + human B2M
mhctransform molecule-to-af3 'HLA-A*02:01' --class I

# 适合 AF3 的 MHC-II alpha/beta 双链
mhctransform molecule-to-af3 \
  'HLA-DRA*01:01/HLA-DRB1*15:01' --class II

# NetMHCpan 兼容的 34-aa descriptor
mhctransform allele-to-pseudo 'HLA-A*02:01' \
  --class I --definition netmhcpan-34-v1
```

反向查询使用 `full-to-allele`、`af3-to-allele` 和 `pseudo-to-allele`。反查会保留全部完全匹配候选，不会把多对一映射伪装成唯一等位基因。详细定义、覆盖范围和审核结论见 [`MHCtransform/README.md`](MHCtransform/README.md) 与 [`MHCtransform/AUDIT_REPORT.md`](MHCtransform/AUDIT_REPORT.md)。

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

`MHCtransform` 的 `af3-soluble-extracellular-v1` 会移除 signal peptide、跨膜区和胞内尾，保留全部胞外折叠域及膜近端连接区；`TCRtransform` 输出成熟的 TCR V+C 胞外结构链。AF3 官方输入格式只要求分别提供各条氨基酸链，并未规定统一裁剪坐标；上述范围是本仓库经过来源核对和测试的工程 profile。

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

本仓库代码采用 [MIT License](LICENSE)。每个子项目的 `NOTICE.md` 和数据目录元数据记录了 IPD-IMGT/HLA、UniProt、tcrdist3、PUFFIN 等上游来源、固定版本和校验值。MIT 许可证只适用于本仓库代码，不会重新许可所捆绑或派生的生物参考数据；再分发数据前请核对各上游来源的适用条款。
