# DNA Restriction Digest Simulator

基于 Python + FastAPI 的 DNA 限制性内切酶平切模拟服务。接收 ACGT 序列、线性/环状标记和酶规则，搜索识别串及其反向互补串（保留重叠匹配），合并所有切点同时切割，返回各酶切点与酶切片段。亦可为多个 DNA 样本设计最小酶组合（电泳面板），使任意两个样本的条带图谱都可区分。`/partial-digest` 在上述匹配/去重规则之上，为每种酶增加独立切开概率 p，精确计算片段长度的期望数量分布与整条分子无切口概率。`/release-fragment` 则反向求解：给定目标区间，从候选酶中挑选能恰好切出该目标片段的最少酶组合。

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

反向设计：给定目标区间，从 1～12 种候选酶中挑选能**恰好切出该目标片段**的最少酶组合。接口地址 `http://localhost:8000/release-fragment`（端口同其他接口，默认 8000，可用 `PORT` 环境变量修改）。

请求体：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sequence` | string | DNA 序列，仅含 A/C/G/T（大小写均可，内部统一大写） |
| `topology` | string | `"linear"`（默认）或 `"circular"` |
| `start` | int | 目标片段起点（零基，≥ 0），沿原序列正向 |
| `end` | int | 目标片段终点（零基，不含），沿原序列正向 |
| `enzymes` | array | 1～12 种候选酶，格式与 `/digest` 相同：`name`（唯一）、`site`（IUPAC）、`cut_offset`（0 ≤ k ≤ 串长，且识别串不得长于序列） |

区间合法性（不满足返回 422）：

- 线性：`0 ≤ start < end ≤ 序列长度`。
- 环状：`start`、`end` 均为 `< 序列长度` 的合法坐标且**不相等**；`end < start` 时目标片段跨首尾（正向绕回原点）。

规则要点：

- 复用 `/digest` 的识别串/反向互补匹配、重叠保留、同酶同点去重与环状切点取模规则。
- 目标两端必须各为**所选酶的切点**或**线性天然末端**（线性 `start = 0` 或 `end = n` 时该端无需切点；环状无天然末端，两端都必须是切点）。
- 所选酶在目标**内部**（两端之间的弧上，不含端点）不得产生任何切点；在目标之外切割不受限制。
- 同一位置被多种酶命中只算**一个**切点（响应中的 `cuts` 为所选酶全部切点去重排序后的列表）。
- 选择满足条件的最少酶组合；酶数相同的方案按**排序后酶名列表的字典序**取首个。
- 线性全长目标（`start = 0` 且 `end = n`）两端均为天然末端，不需要任何酶：返回空酶列表。
- 无解时返回明确的无解结果：`feasible = false`，`enzymes`/`cuts` 为空，`fragment` 为 `null`（HTTP 状态仍为 200）。

响应字段：

| 字段 | 说明 |
|---|---|
| `feasible` | 是否存在满足条件的酶组合 |
| `enzymes` | 所选酶名（排序后）；无解时为空 |
| `fragment` | 释放出的目标片段（`start`、`length`、原序列方向的 `sequence`）；无解时为 `null` |
| `cuts` | 所选酶在分子上的实际切点（去重、升序，含目标之外的切点）；无解时为空 |

示例（线性，单酶即可切出目标片段）：

```bash
curl -s -X POST http://localhost:8000/release-fragment \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "ACGTGAATTCGGGAATTCAA",
        "topology": "linear",
        "start": 5,
        "end": 9,
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}
        ]
      }'
```

响应（EcoRI 在 5 和 9 的切点正好为目标两端，13、17 为目标之外的实际切点）：

```json
{
  "feasible": true,
  "enzymes": ["EcoRI"],
  "fragment": {"start": 5, "length": 4, "sequence": "AATT"},
  "cuts": [5, 9, 13, 17]
}
```

环状跨首尾示例（目标从 5 正向绕回 1）：

```bash
curl -s -X POST http://localhost:8000/release-fragment \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "GAATTC",
        "topology": "circular",
        "start": 5,
        "end": 1,
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}
        ]
      }'
```

```json
{
  "feasible": true,
  "enzymes": ["EcoRI"],
  "fragment": {"start": 5, "length": 2, "sequence": "CG"},
  "cuts": [1, 5]
}
```

无解示例（目标两端都不是切点）：

```bash
curl -s -X POST http://localhost:8000/release-fragment \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "ACGTGAATTCGGGAATTCAA",
        "start": 6,
        "end": 8,
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}
        ]
      }'
```

```json
{
  "feasible": false,
  "enzymes": [],
  "fragment": null,
  "cuts": []
}
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
  release.py   # 片段释放设计：目标区间内部切点过滤、必需端点覆盖、最少酶组合枚举
  schemas.py   # 请求/响应模型与校验（含 /partial-digest 的 p 与 2000 bp 限制、/release-fragment 的区间校验）
tests/
  test_digest.py
  test_partial_digest.py
  test_panel.py
  test_release.py
Dockerfile
docker-compose.yml
requirements.txt
```
