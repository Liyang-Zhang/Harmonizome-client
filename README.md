# Harmonizome-client

本仓库包含两个主要用途：
- 通过 CLI 脚本下载 Harmonizome 的公开数据集（REST/KG API、批量数据）。
- 使用本地下载的 edge 文件对基因列表进行离线注释（每个输入基因一行，追加注释列）。

## 目录结构
- `harmonizome_cli.py`：下载 Harmonizome 数据集文件（`-d` 数据集 key，`-t` 文件类型）。
- `harmonizome_api_cli.py` / `harmonizome_api_update.py`：调用 Harmonizome REST + KG API 的命令行工具（更新版含 export-genes）。
- `config/`
  - `databasesheet.csv`：手工配置的数据库名称与 edge 文件 URL。
  - `local_datasets.yml`：注释用的本地数据集配置（当前包含 GTEx/HPA/GO BP/ClinVar/OMIM 等 edge 文件，可自行增删）。
- `downloaded_edges/`：按配置从网页下载的 edge 文件存放处（每个数据集一个子目录）。
- `harmonizome_data/`：通过 `harmonizome_cli.py` 下载的文件存放处（之前用 API/CLI 拉取的版本）。
- `annotate_local.py`：离线注释脚本（读取基因列表 + 本地 edge 文件，输出 TSV）。
- 其他：`Harmonizome.pdf` 文档、`sample.jsonl` 示例、`annotations.tsv` 最新注释结果等。

## 环境依赖
- Python 3.8+
- `requests`（已用于下载/HTTP）
- `PyYAML`（读取 YAML 配置；若缺失可以改用 JSON 配置）

安装示例：
```bash
python3 -m pip install requests PyYAML
```

## 下载数据
- 根据 `config/databasesheet.csv` 的 URL，示例下载（已放在 `downloaded_edges/`）：
  - GTEx Tissue Gene Expression Profiles
  - HPA Tissue Gene Expression Profiles
  - ClinVar Gene-Phenotype Associations 2025
  - OMIM Gene-Disease Associations
  - GO Biological Process Annotations 2025
- 如需新增数据集，可在 `databasesheet.csv` 添加名称和 edge 文件 URL，自行下载或编写脚本读取。
- 也可用 `harmonizome_cli.py` 通过数据集 key 下载，例如：
  ```bash
  python3 harmonizome_cli.py download -d gtextissue -t gene_attribute_edges.txt.gz -o harmonizome_data
  ```

## 离线注释（annotate_local.py）
- 输入：基因列表
  - 支持命令行 `--genes` 直接列出基因，或
  - 通过 `--genes-file` 指定一个 TSV 文件，文件需包含列名为 `gene` 或 `基因`，其余列会原样保留。如列名不规范，可用 `--gene-col-index` 指定基因列序号（1-based，优先）。
- 配置：`config/local_datasets.yml` 定义每个数据集的本地路径、基因列名、输出模式。
  - 默认每个数据集输出两列：`<name>_json`（匹配记录 JSON 序列化）和 `<name>_count`（命中条数），防止列数爆炸。
  - 支持 `mode: tissues`（如当前 GTEx/HPA Protein），按配置的组织列表输出多列，列名形如 `<name>_<tissue>`，值为 1/-1/占位符（不输出 count 列）。
  - 支持 `mode: cells`（当前 ENCODE_HM），按配置的 cell × mark 组合输出列，列名形如 `<name>_<cell>_<mark>`，值为 1/-1/占位符（可选 marks 过滤）。
- 输出：TSV，行数与输入行数一致，原列 + 注释列。

示例：
```bash
# 使用命令行基因
python3 annotate_local.py --genes EGFR TP53 KRAS ERBB2 \
  --config config/local_datasets.yml \
  --out annotations.tsv

# 使用包含“基因”列的 TSV 文件
python3 annotate_local.py \
  --genes-file ./D121X251020ECS0001AR001_1_GATK_ACMG_exon20bp_panel.txt \
  --config config/local_datasets.yml \
  --out annotations.tsv
# 如输入 TSV 的基因列不是常规名称，也可指定列号（1-based）：
# --gene-col-index 3
```

## 扩展/定制
- 新增数据集：在 `config/local_datasets.yml` 里添加条目，指定 `name/path/gene_field/fields/output`。
- 输出模式：默认 `json`；如需扁平化/拼接可将 `output.mode` 设为 `join` 并配置 `join_field/sep`，或自行扩展脚本的聚合逻辑。
- 如果想对未阈值化或标准化分值的文件进行注释，只需在配置中指向对应文件，并在 `fields` 中选择要保留的字段。
- 当前示例配置包含：GTEx_Tissue（mode=tissues）、HPA_Protein（mode=tissues）、ENCODE_HM（mode=cells，细胞/标记列表在 YAML 有中文注释），不输出 count 列；如需其他库可在 YAML 中添加/恢复。

## 已知差异
- 网页端下载的部分 edge 文件包含标准化分值且命名更细粒度；API/CLI 下载的版本多为阈值化 ±1。根据需求选择对应文件并在配置中调整路径。
