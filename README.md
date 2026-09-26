# DNA Restriction Digest Simulator

基于 Python + FastAPI 的 DNA 限制性内切酶平切模拟服务。接收 ACGT 序列、线性/环状标记和酶规则，搜索识别串及其反向互补串（保留重叠匹配），合并所有切点同时切割，返回各酶切点与酶切片段。亦可为多个 DNA 样本设计最小酶组合（电泳面板），使任意两个样本的条带图谱都可区分。`/partial-digest` 在上述匹配/去重规则之上，为每种酶增加独立切开概率 p，精确计算片段长度的期望数量分布与整条分子无切口概率。`/release-fragment` 为反向设计：给定目标片段区间，从候选酶中选出能恰好切出该片段的最少酶组合。

## 启动

单条命令（需要 Docker 与 Docker Compose）：

```bash
docker compose up --build
```

镜像构建时安装依赖，容器启动时会再次执行 `pip install -r requirements.txt` 确认依赖就绪，然后启动 uvicorn。

## 配置

| 环境变量 | 默认值 | 说明 |
|---|---|---|
| `PORT` | `8000` | 宿主机映射端口，如 `PORT=9000 docker compose up --build` |

## API 地址

- 服务地址： `http://localhost:8000`
- 交互式文档（Swagger UI）： `http://localhost:8000/docs`
- 健康检查： `GET http://localhost:8000/health` → `{"status": "ok"}`
- 酶切模拟： `POST http://localhost:8000/digest`
- 概率酶切模拟： `POST http://localhost:8000/partial-digest`
- 酶面板设计： `POST http://localhost:8000/design-panel`
- 目标片段释放设计： `POST http://localhost:8000/release-fragment`

## 接口说明

### `POST /digest`

请求体：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sequence` | string | DNA 序列，仅含 A/C/G/T（大小写均可，内部统一大写） |
| `topology` | string | `"linear"`（默认）或 `"circular"` |
| `enzymes` | array | 酶规则列表，每条含 `name`（唯一）、`site`（IUPAC 识别串）、`cut_offset`（串内切割偏移 k，0 ≤ k ≤ 串长） |

规则要点：

- 同时搜索识别串及其反向互补串，保留重叠匹配；正向切点 = 起点 + k，反向切点 = 起点 + 串长 − k。
- 环状序列允许跨首尾的匹配，切点按序列长度取模；识别串不得长于序列（否则返回 422）。
- 切点为零基坐标，同酶同点去重；所有酶的切点合并后同时切割。
- 片段按起点排序，含 `start`、`length` 和原序列方向的碱基串；环状末段绕回首段。
- 无切点时返回原序列作为唯一片段；不输出空片段。

示例（线性，双酶切）：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "AGGATCCGAATTCA",
        "topology": "linear",
        "enzymes": [
          {"name": "BamHI", "site": "GGATCC", "cut_offset": 1},
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}
        ]
      }'
```

响应：

```json
{
  "cuts": [
    {"name": "BamHI", "cuts": [2, 6]},
    {"name": "EcoRI", "cuts": [8, 12]}
  ],
  "fragments": [
    {"start": 0,  "length": 2, "sequence": "AG"},
    {"start": 2,  "length": 4, "sequence": "GATC"},
    {"start": 6,  "length": 2, "sequence": "CG"},
    {"start": 8,  "length": 4, "sequence": "AATT"},
    {"start": 12, "length": 2, "sequence": "CA"}
  ]
}
```

环状示例（末段绕回首段）：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "GAATTC",
        "topology": "circular",
        "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}]
      }'
```

```json
{
  "cuts": [{"name": "EcoRI", "cuts": [1, 5]}],
  "fragments": [
    {"start": 1, "length": 4, "sequence": "AATT"},
    {"start": 5, "length": 2, "sequence": "CG"}
  ]
}
```

### `POST /partial-digest`

概率酶切：每种酶在自己的**每个**切点处以概率 `p`（0～1）独立切开；同一位置被多种酶命中时，任一种切开即形成切口。服务**精确计算**（解析概率乘积，不采样、不枚举切口子集），以有效切点数 K 计，分布计算时间不超过 O(K²)。

请求体（与 `/digest` 相同，另加 `p`，序列最长 **2000 bp**）：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sequence` | string | DNA 序列，仅含 A/C/G/T，长度 1～2000（超长返回 422） |
| `topology` | string | `"linear"`（默认）或 `"circular"` |
| `enzymes` | array | 酶规则列表，每条含 `name`（唯一）、`site`、`cut_offset` 以及 `p`（0～1，越界返回 422）；空列表视为不切 |

规则要点：

- 切点匹配、反向互补、重叠保留、同酶同点去重、环状切点取模等规则与 `/digest` 完全一致。
- 各酶在各位置的切开事件相互独立；同一物理位置的多酶切开概率合并为 `1 − Π(1 − p)`。
- **线性分子两端（坐标 0 与 n）不是随机切口**；只有内部坐标可以形成切口。
- `p = 0` 的酶不参与切割；`p = 1` 时结果即确定性酶切，可与 `/digest` 互验。
- 返回按长度升序的**片段期望数量**（同长度累加，计数可为小数；零值省略）以及**整条分子无切口概率**。
- 环状分子：无切口和恰有一个切口都产生一个全长片段（期望相加），但 `no_cut_probability` 只计无切口情形。

示例（环状，p = 0.3）：

```bash
curl -s -X POST http://localhost:8000/partial-digest \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "GAATTCGAATTC",
        "topology": "circular",
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 0.3}
        ]
      }'
```

响应：

```json
{
  "expected_fragments": [
    {"length": 2,  "count": 0.18},
    {"length": 4,  "count": 0.18},
    {"length": 6,  "count": 0.252},
    {"length": 8,  "count": 0.0882},
    {"length": 10, "count": 0.0882},
    {"length": 12, "count": 0.6517}
  ],
  "no_cut_probability": 0.2401
}
```

其中无切口概率 = 0.7⁴ = 0.2401；全长（12 bp）期望 0.6517 = 无切口 0.2401 + 恰有一个切口 4 × 0.3 × 0.7³ = 0.4116。

### `POST /design-panel`

为 2～12 个 DNA 样本挑选能区分所有样本对的最少酶组合（电泳面板）。

请求体：

| 字段 | 类型 | 说明 |
|---|---|---|
| `samples` | array | 2～12 个样本，每个含 `name`（唯一）、`sequence`（仅 A/C/G/T，大小写均可）、`topology`（`"linear"` 默认 / `"circular"`） |
| `enzymes` | array | 1～12 种候选酶，格式与 `/digest` 相同：`name`（唯一）、`site`（IUPAC）、`cut_offset`（0 ≤ k ≤ 串长，且识别串不得长于任一样本序列） |
| `tolerance_bp` | int | 条带容差（bp），非负整数 |

规则要点：

- 复用 `/digest` 的酶切逻辑：每个样本分别用每种候选酶做单酶切，每种酶独立成一条泳道。
- 泳道条带 = 片段长度升序列表，重复条带保留。
- 两条泳道条带数量相同且对应长度差均 ≤ `tolerance_bp` 时视为不可区分，否则该酶可区分这对样本。
- 选择能区分全部样本对的最少酶组；同数量方案按排序后的酶名列表字典序取最小。
- 无解（全部候选酶合起来仍无法区分某些样本对）时：`enzymes`、`lanes`、`pairs` 为空，`undistinguishable_pairs` 列出这些样本对。

响应字段：

| 字段 | 说明 |
|---|---|
| `enzymes` | 所选酶名（排序后）；无解时为空 |
| `lanes` | 每个样本在所选泳道中的条带（`bands`：酶名 → 升序长度列表），按样本输入顺序 |
| `pairs` | 每对样本及可区分它们的所选酶 |
| `undistinguishable_pairs` | 全部候选酶仍无法区分的样本对（有解时为空） |

示例：

```bash
curl -s -X POST http://localhost:8000/design-panel \
  -H 'Content-Type: application/json' \
  -d '{
        "samples": [
          {"name": "wild-type", "sequence": "AAGAATTCTT"},
          {"name": "mutant",    "sequence": "AAGGATCCTT"}
        ],
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
          {"name": "BamHI", "site": "GGATCC", "cut_offset": 1}
        ],
        "tolerance_bp": 0
      }'
```

响应（两种酶单独都能区分该样本对，按字典序取 BamHI）：

```json
{
  "enzymes": ["BamHI"],
  "lanes": [
    {"sample": "wild-type", "bands": {"BamHI": [10]}},
    {"sample": "mutant",    "bands": {"BamHI": [3, 3, 4]}}
  ],
  "pairs": [
    {"samples": ["wild-type", "mutant"], "enzymes": ["BamHI"]}
  ],
  "undistinguishable_pairs": []
}
```

### `POST /release-fragment`

反向设计：给定目标片段区间，从至多 12 种候选酶中找出能切出**恰好该片段**的最少酶组合。

请求体：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sequence` | string | DNA 序列，仅含 A/C/G/T（大小写均可，内部统一大写） |
| `topology` | string | `"linear"`（默认）或 `"circular"` |
| `start` | int | 目标片段起点（零基坐标，≥ 0），沿原序列正向 |
| `end` | int | 目标片段终点（零基坐标，≥ 0） |
| `enzymes` | array | 0～12 种候选酶，格式与 `/digest` 相同：`name`（唯一）、`site`（IUPAC）、`cut_offset`（识别串不得长于序列） |

规则要点：

- 目标片段沿原序列正向从 `start` 到 `end`。线性要求 `0 ≤ start < end ≤ 序列长度`；环状要求 `start`、`end` 为 `[0, 序列长度)` 内两个**不同**坐标，允许跨首尾（`start > end` 时绕回原点）。区间不合法返回 422。
- 识别位点搜索、反向互补、重叠保留、切点偏移与环状取模规则与 `/digest` 完全一致。
- 释放条件：目标两端都是所选酶产生的切点（线性分子的天然末端 0 与 n 也算），且所选酶在目标**内部**不产生任何切点；目标之外的切点不受限制。
- 同一点位被多种酶切割只算一个切点；`cuts` 为所选酶合并去重后的全部实际切点（升序）。
- 按子集大小递增枚举，第一个可行组合即为答案；等量方案按排序后的酶名列表字典序取首个。线性全长目标（`start=0`、`end=序列长度`）不需要任何酶，返回空酶列表。
- 无解时仍返回 HTTP 200：`solved` 为 `false`，`enzymes`/`cuts` 为空、`fragment` 为 `null`。

响应字段：

| 字段 | 说明 |
|---|---|
| `solved` | 是否找到可行酶组合 |
| `enzymes` | 所选酶名（排序后）；无解时为空 |
| `fragment` | 目标片段（`start`、`length`、`sequence`，环状可跨首尾）；无解时为 `null` |
| `cuts` | 所选酶合并去重后的实际切点（升序）；无解时为空 |

示例（线性，单酶即可切出）：

```bash
curl -s -X POST http://localhost:8000/release-fragment \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "AAGAATTCTTGGATCCAA",
        "topology": "linear",
        "start": 3,
        "end": 7,
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
          {"name": "BamHI", "site": "GGATCC", "cut_offset": 1}
        ]
      }'
```

响应（BamHI 的切点 11/15 不在目标两端，单选 EcoRI 即可）：

```json
{
  "solved": true,
  "enzymes": ["EcoRI"],
  "fragment": {"start": 3, "length": 4, "sequence": "AATT"},
  "cuts": [3, 7]
}
```

环状跨首尾示例：

```bash
curl -s -X POST http://localhost:8000/release-fragment \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "GAATTC",
        "topology": "circular",
        "start": 5,
        "end": 1,
        "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}]
      }'
```

```json
{
  "solved": true,
  "enzymes": ["EcoRI"],
  "fragment": {"start": 5, "length": 2, "sequence": "CG"},
  "cuts": [1, 5]
}
```

无解示例（EcoRI 在目标内部 9、13 处还有切点）：

```bash
curl -s -X POST http://localhost:8000/release-fragment \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "ACGTGAATTCGGGAATTCAA",
        "start": 5,
        "end": 17,
        "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}]
      }'
```

```json
{"solved": false, "enzymes": [], "fragment": null, "cuts": []}
```

## 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试（需额外安装 `pytest` 与 `httpx`）：

```bash
python3 -m pytest tests/ -q
```

## 项目结构

```
app/
  main.py      # FastAPI 应用与路由（/health、/digest、/partial-digest、/design-panel、/release-fragment）
  digest.py    # 核心酶切逻辑：IUPAC 匹配、切点计算、片段生成、概率酶切期望分布（O(K²)）
  panel.py     # 面板设计：单酶切条带、可区分性判定、最小酶组选择
  release.py   # 片段释放设计：目标区间判定、最少酶组合搜索（按大小递增 + 字典序破平）
  schemas.py   # 请求/响应模型与校验（含 /partial-digest 的 p 与 2000 bp 限制、/release-fragment 的区间与 12 酶上限）
tests/
  test_digest.py
  test_partial_digest.py
  test_panel.py
  test_release.py
Dockerfile
docker-compose.yml
requirements.txt
```
