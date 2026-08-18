# bioinfo05 — OA 与 OP 共病的分区特异转录组收敛：RNA 剪接程序轴

> 多组学、无偏的骨关节炎（OA）–骨质疏松（OP）共病发现项目。
> 核心论点：**OA 与 OP 的转录组收敛是"分区特异"的**——软骨（间充质）与骨髓基质细胞（间充质）共享一条一致下调的 RNA 剪接程序，而循环单核（造血）分区走向相反，因此"按疾病标签合并"的分析会把信号抵消为阴性。

---

## 1. 核心发现

| # | 层次 | 结论 | 关键证据 |
|---|------|------|----------|
| 1 | 基因 | OA 与 OP **无**共享差异表达基因 | 仅 10 个交集，超几何 P = 0.67；meta z 全基因组相关 −0.020 (P = 0.16) |
| 2 | 通路（合并） | 按疾病标签 pooling 时**无**收敛 | 富集 0.62，置换 P = 0.93 |
| 3 | 通路（分区） | 间充质分区**强收敛** | OA 软骨 + OP BM-MSC 富集 3.33 (P = 0.030)；家族级 66 vs 55.1 (P = 0.004) |
| 4 | 通路（对照） | 造血分区反而**低于**零分布 | OA 软骨 vs OP 单核 P = 0.997，rho = −0.065 |
| 5 | 基因轴 | **49 个核心剪接因子**在间充质两组织一致下调 | 含 SF3B1、HNRNPU、SRRM1、PRPF8、CDC5L；仅 45% 在 OP 单核也下调（binom P = 0.57）→ 分区特异 |
| 6 | 单细胞 | 剪接位移**未定位于特定软骨亚型**（阴性，欠力） | 3 vs 3 供体，pseudobulk 后所有亚型 P > 0.45 |
| 7 | 机器学习 | 跨分区分类器**不迁移**（阴性，欠力） | AUC 0.30 (n = 9)，零分布 sd = 0.245 |
| 8 | 遗传学 | GWAS Catalog 提供**观察性**收敛证据 | 70 候选基因中 17 个基因座在骨/关节性状显著（114 条记录），最强为 heel BMD（CDC5L 位点 P < 1e-100） |

**GO 功能注释（正式，以 10,373 共同基因为背景）**：WGCNA 跨病种 21 共享基因**无**显著 BP 通路（FDR > 0.05），仅作描述性；剪接程序 205 基因注释到 95 个 GO BP（RNA splicing FDR = 3.1e-198 等）。

---

## 2. 数据来源

公共 GEO 数据集（全部去标识化）：

- **OA 软骨（4）**：GSE114007、GSE57218、GSE117999、GSE169077
- **OP 循环单核（2）**：GSE56815、GSE7158
- **OP BM-MSC（1，9 样本）**：GSE35958
- **单细胞软骨图谱**：GSE324993（37,453 细胞，3 健康 + 3 OA 供体）
- **遗传**：GWAS Catalog alt-full associations release（查询骨/关节性状）

7 个 bulk 队列经批次校正后共同量化 **10,373** 个基因，作为所有下游分析的 universe。

> **数据文件说明**：本仓库随附分析所需的处理后数据（`01_data/processed/`，表达矩阵与表型）与基因集（`01_data/genesets/`，MSigDB v7.5.1）。**原始 GEO 系列矩阵、单细胞原始数据、芯片注释文件及 GWAS Catalog 全量关联表未随仓库分发**（体积过大），可从公共源按上述 accession / release 自行下载后置于对应目录运行脚本。

### `01_data/gwas/` 使用须知

| 文件 | 状态 | 说明 |
|------|------|------|
| `gwas_assoc_full.zip`（70 MB） | ✅ 完整可用 | GWAS Catalog 官方全量关联表（alt-full associations release），需自行从 https://www.ebi.ac.uk/gwas/docs/file-downloads 下载后置于 `01_data/gwas/`。运行 `17_genetics_gwas_catalog.py` 会自动解压（解压后 701 MB TSV 为可再生中间件，不入库） |
| `eBMD_Morris2019_GCST006979.txt.gz`<br>`Fracture_Morris2019_GCST006980.txt.gz`<br>`OA_Tachmazidou2019_GCST007093.txt.gz` | ⚠️ **下载不完整** | 各约 15 MB，为 MR 尝试期间中断的部分文件。**请勿用于任何分析**——截断的 sumstats 会导致无效的工具变量选择。保留仅为记录 MR 未能完成的过程

---

## 3. 分析流程（`02_code/`）

| 脚本 | 内容 |
|------|------|
| `01_preprocess_OP.R` | OP 队列预处理 |
| `02_harmonise_OA_and_QC.R` | OA 协调与质控 |
| `03_DEG.R` | 差异表达（limma / limma-voom） |
| `04_SON_diagnostic.R` | SON 候选 hub 诊断 |
| `05_*` | 无偏发现：WGCNA (`05a`)、GSVA 收敛 (`05b`)、收敛对照 (`05c`)、标签置换零模型 (`05d`)、分区收敛 (`05e`)、剪接轴 (`05f`)、置换/循环性 (`05g`)、通路家族检验 (`05h`) |
| `06_singlecell.R` | 单细胞分析（Seurat） |
| `07_ml_compartment_transfer.R` | 跨分区 ML 迁移 |
| `08_verify_references.py` / `08b_*` | 参考文献真实性核查（Europe PMC / Crossref） |
| `09_figures_1_3.R` / `10_figures_4_6.R` | Figure 1–6 |
| `11_tables_word.R` | 表格生成 |
| `15_functional_annotation.R` | GO 功能注释 |
| `16_write_paper_final.R` | **终稿论文**（两遍渲染 [n] 引用，flextable 直接写表） |
| `17_genetics_gwas_catalog.py` | **GWAS Catalog 遗传支持**（自动从 zip 解压后流式检索 701 MB TSV） |
| `18_finalize_docx.py` | 修复 flextable 缺 `w:tblGrid` 的缺陷，保证文档可被 Word 打开 |
| `_validate_project.R` / `19_validate_project.py` | 交付物完整性校验 |
| `20_compartment_robustness.py` | **双 OP 单核队列区室稳健性检验**（GSE56815+GSE7158）：协调剪接程序在两单核队列均不存在；SON 跨 7 队列实测 |
| `21_validate_v2.py` | 针对作者/基金/声明/Table 3 符号/Table 2a "vs"/Abstract FDR/Figure 编号的专项校验 |
| `22_make_cover_letter.py` | **生成 Cover_Letter.docx**（含显式 "re-analysis resolves contradictions in Chen 2024 / Wang 2025 / Jiang 2026" 句，预印/投稿用） |
| `23_validate_v3.py` | 10 项小修专项校验（Table 3 负号 / Table 2a vs / Table 4 脚注 11·49 / Figure 6 / Abstract FDR / 引用[3][8] / 图编号 / Methods 公开仓库） |
| `24_validate_v4.py` | v4 校验：作者上标 ¹¹¹²\* + 17 项回归（含 Ref [3] 干净核对，DOI 正确、无 PMID 内联） |

> 注：`12_gwas_catalog_top_genes.py`、`13_write_paper.R`、`14_combine_paper.R` 为早期版本，因假阳性文献/altChunk 损坏已弃用，保留作中间产物。

---

## 4. 主要交付物

- **`06_paper/Manuscript_final_v4.docx`** — **当前最新发表级终稿**（6 图 + 5 表，26 篇已核实参考文献，两遍渲染编号；作者上标 ¹¹¹²\* 已修正，单位 1/2 分离；含作者/基金/声明、Table 3 符号、11/49 脚注统一、公开仓库、re-analysis Cover Letter）。若被外部预览锁住无法覆盖，直接打开此文件即可。历史版本：`Manuscript_final_v2.docx`（含作者/基金/声明）、`Manuscript_final_ready.docx`
- **`04_figures/Fig1–6.{png,tiff}`** — 发表级图表（tiff 为印刷格式）
- **`05_tables/`** — Table 1（队列）、2a/2b（收敛与家族级）、3（剪接通路）、49 基因表
- **`03_results/`** — `deg/`、`gsva/`、`wgcna/`、`functional/`（GO）、`genetics/`（GWAS Catalog）、`ml/`、`sc/` 等
- **`03_results/functional/compartment_robustness.txt`** / **`SON_per_cohort.csv`** — 双单核队列区室稳健性结果（20_compartment_robustness.py 输出）
- **`06_paper/references_verified.csv`** — 54 条已核实文献（26 条被引用）
- **`00_docs/00_original_research_idea.docx`** — 原始研究思路文档（溯源存档）

---

## 5. 复现

```bash
cd D:/bioinfo05
# R 4.6.1
"path/to/Rscript.exe" 02_code/16_write_paper_final.R
# Python 3.13 (venv) — 修复 docx 表格网格，保证可打开
python 02_code/18_finalize_docx.py
# 可选：遗传支持 / 校验
python 02_code/17_genetics_gwas_catalog.py
python 02_code/19_validate_project.py
```

随机种子固定；bulk 分析用 R 4.6.1，docx 修复与遗传检索用 Python 3.13 venv。

---

## 6. 方法学严格性

- **三套零分布**：通路标签置换（2,000）、病例标签置换（500）、基因集家族 Jaccard 聚类（2,168 家族）消除伪独立。
- **阳性对照**：两个不同平台的 OA 队列强收敛（富集 3.06, P < 0.001），验证流程敏感度。
- **正式 GO**：以 10,373 共同基因为背景，避免全基因组虚高。
- **文献零虚构**：48 条参考文献经 Europe PMC / Crossref 逐条核实（含 Benjamini–Hochberg 1995 原始方法学文献）。

---

## 7. 已知局限（已在文内诚实报告）

1. 仅 1 个 BM-MSC 队列（9 样本），间充质轴不对称地依赖小样本。
2. 单细胞 3 vs 3 供体，检验力不足——阴性结果约束解释而非推翻主发现。
3. **未完成 Mendelian randomisation**：全量 GWAS summary statistics 下载速率过低（~45 KB/s），无法支持有效工具变量选择；以 GWAS Catalog 观察性证据替代，方向性因果未定。
4. 剪接因子转录丰度为**间接** proxy，未实测 splice-junction 水平 isoform 使用。

---

## 8. 作者与基金（已填写）

论文末 Declarations / Funding / Author contributions 已填写完成：

- **作者**：Wei Yuwei¹, Li Peng¹, Chen Kai¹, Zhang Dong²,*（¹ 安徽省淮北市濉溪县中医院骨伤科；² 江苏省人民医院重庆医院骨科中心）
- **通讯作者**：Zhang Dong, MD, Chief Physician. Email: zhangdht@126.com. ORCID: 0009-0006-9000-273X
- **基金**：淮北市科技局自筹科技计划项目（No. 2025HK044）
- 本仓库不含未脱敏的个人信息。
