# **KCNQ1 F335A 突变对 ML277 结合能影响的 FEP 计算工作流**

### **目标**
使用自由能微扰（FEP）方法，计算 KCNQ1 蛋白的 F335A 点突变对配体 ML277 结合自由能的影响（ΔΔG）。

### **核心策略：热力学循环**
我们将通过计算两个独立的、计算上可行的“炼金术”突变路径的自由能（ΔG）来达成目标：
1.  **ΔG_apo**: 在无配体的 Apo 蛋白中，将 F335 突变为 A335。
2.  **ΔG_complex**: 在与 ML277 结合的 Holo 复合物中，进行同样的 F335A 突变。

最终的结合自由能变化为： **\[ ΔΔG_{binding} = ΔG_{complex} - ΔG_{apo} \]**
---------------
您好！这是一个非常有趣且有价值的计算化学项目。使用 FEP (Free Energy Perturbation) 来计算点突变对配体结合自由能的影响是该领域的黄金标准方法之一。

根据您的描述和已有文件，您已经为野生型（WT）的 `apo` 状态蛋白打下了良好的基础。接下来，我们将利用热力学循环（Thermodynamic Cycle）来高效地计算结合自由能的变化（ΔΔG）。

### 核心策略：热力学循环

我们不直接计算 ML277 与两种蛋白（野生型和突变型）的结合能（ΔG_binding_WT, ΔG_binding_MUT），因为模拟配体的结合/解离过程非常耗时。取而代之，我们计算两个非物理过程（炼金术突变）的自由能变化：

1. **Apo 状态下的突变 (Leg 1):** 在没有配体的情况下，将 KCNQ1 的 335 位残基从苯丙氨酸（F）突变为丙氨酸（A）。这个过程的自由能变化是 **ΔG_apo**。
  `KCNQ1_WT (apo) → KCNQ1_F335A (apo)`
  
2. **Holo 状态下的突变 (Leg 2):** 在与配体 ML277 结合的状态下，进行同样的 F→A 突变。这个过程的自由能变化是 **ΔG_complex**。
  `KCNQ1_WT-ML277 (complex) → KCNQ1_F335A-ML277 (complex)`
  

根据热力学循环，我们最终关心的结合自由能变化可以通过以下公式得到：

**\[ ΔΔG_{binding} = ΔG_{binding}^{MUT} - ΔG_{binding}^{WT} = ΔG_{complex} - ΔG_{apo} \]**

这个差值（ΔΔG）直接反映了 F335A 突变对 ML277 结合亲和力的影响。如果 ΔΔG < 0，说明突变增强了结合；如果 ΔΔG > 0，说明突变削弱了结合。

---

### 重要前提：四聚体 vs. 单体

KCNQ1 是四聚体，有四个结合口袋。在单个 Colab 环境中模拟整个四聚体（超过20万个原子）的 FEP 计算量极其庞大，几乎不可能完成。

**强烈建议**：从一个单体（protomer）开始模拟。由于四个亚基是相同的，并且如果结合口袋之间没有直接的变构效应，单个口袋的突变可以很好地代表整体情况。您需要从您的四聚体结构中提取出一个包含完整结合口袋的单体以及与之结合的单个 ML277 分子。

- **Apo 状态**: KCNQ1 单体
- **Holo 状态**: KCNQ1 单体 + 1个 ML277

后续的计划将基于**模拟单个亚基**的前提。

---

### 下一步模拟计划

您已经完成了野生型 Apo 蛋白的 MD，这是一个很好的开始。现在我们需要为 FEP 计算准备和运行两个核心“分支”：Apo 突变和 Complex 突变。

下表总结了后续的步骤：

| 阶段/步骤 | 任务描述 | 主要工具 | 预计耗时 (Colab GPU) | 备注/前置条件 |
| --- | --- | --- | --- | --- |
| **1. 系统准备** | 1. 从四聚体中提取**单体**结构 (WT-apo, WT-complex)。<br>2. 使用 `pmx mutate` 生成 F335A 突变的混合拓扑结构。<br>3. 为配体 ML277 生成 GROMACS 拓扑和参数 (GAFF力场)。 | `GROMACS`, `pmx`, `acpype`/`Antechamber` | 1-2 小时 | 拥有 WT-apo 和 WT-complex 的单体 PDB 文件。 |
| **2. Apo-FEP (Leg 1)** | **对 Apo 蛋白进行 F→A 突变模拟**<br>1. 对 `pmx` 生成的每个 λ 窗口进行能量最小化 (EM)。<br>2. 对每个 λ 窗口进行 NVT 平衡。<br>3. 对每个 λ 窗口进行 NPT 平衡。<br>4. 对每个 λ 窗口运行独立的 Production MD (FEP 采样)。 | `pmx`, `pyautoFEP`, `GROMACS` | **2-4 天/分支** | 已完成系统准备。这是计算的主要瓶颈。 |
| **3. Complex-FEP (Leg 2)** | **对 Holo (蛋白-配体复合物) 进行 F→A 突变模拟**<br>1. 对 `pmx` 生成的每个 λ 窗口进行 EM。<br>2. 对每个 λ 窗口进行 NVT 平衡。<br>3. 对每个 λ 窗口进行 NPT 平衡。<br>4. 对每个 λ 窗口运行独立的 Production MD (FEP 采样)。 | `pmx`, `pyautoFEP`, `GROMACS` | **2-4 天/分支** | 已完成系统准备。可与 Apo-FEP 并行运行。 |
| **4. 数据分析** | 1. 收集所有 λ 窗口的 `dhdl.xvg` 文件。<br>2. 使用 `pmx analyse` (或 `gmx bar`) 计算 ΔG_apo 和 ΔG_complex。<br>3. 计算最终的 ΔΔG_binding = ΔG_complex - ΔG_apo。 | `pmx`, `gmx bar` | < 1 小时 | 已完成 Apo 和 Complex 的 FEP 模拟。 |

---

### 工具和代码实现 (`pmx` & `pyautoFEP`)

`pmx` 是核心的突变引擎，而 `pyautoFEP` 是一个非常方便的自动化脚本，可以为您管理和运行 FEP 计算中繁琐的每一步。

**步骤 1: 生成突变 (使用 `pmx mutate`)**

首先，您需要一个包含突变信息的 `mutate.dat` 文件。内容很简单：

```
F 335 A
```

然后，在您的 Colab 环境中，针对 Apo 和 Complex 状态分别运行：

```python
# 假设您的野生型单体蛋白文件为 KCNQ1_WT_mono.pdb
!pmx mutate -f KCNQ1_WT_mono.pdb -o KCNQ1_F335A_apo -ff amber99sb-ildn -water tip3p -dat mutate.dat

# 假设您的野生型复合物文件为 KCNQ1_ML277_WT_mono.pdb
!pmx mutate -f KCNQ1_ML277_WT_mono.pdb -o KCNQ1_F335A_complex -ff amber99sb-ildn -water tip3p -dat mutate.dat
```

这会为 Apo 和 Complex 两个分支分别创建一个文件夹 (`KCNQ1_F335A_apo` 和 `KCNQ1_F335A_complex`)，里面包含了所有 λ 窗口（通常是 12-24 个）的初始结构和拓扑。

**步骤 2 & 3: 运行 FEP 模拟 (使用 `pyautoFEP`)**

`pyautoFEP` 可以自动化处理 EM、平衡和 Production MD。您需要为它准备一个输入文件（例如 `fep.in`），内容大致如下：

```ini
[fep]
; System information
fep_folder = KCNQ1_F335A_apo/
gromacs_path = gmx

; Simulation parameters
em_nsteps = 5000
nvt_nsteps = 50000  ; 100 ps
npt_nsteps = 50000  ; 100 ps
md_nsteps = 1000000 ; 2 ns per window

; Execution parameters
n_threads = 2
gpu_id = 0
```

然后运行它：

```python
# 运行 Apo 分支
!pyautoFEP -i fep_apo.in

# 运行 Complex 分支 (修改 fep.in 中的 fep_folder)
!pyautoFEP -i fep_complex.in
```

**注意：Colab 限制**
Colab 的运行时长有限（免费版约12小时）。FEP 计算通常需要几天。您必须：

1. **频繁使用检查点（Checkpointing）**：`pyautoFEP` 和 `gmx mdrun` 都支持从检查点恢复。您需要分段运行模拟。
2. **管理数据**：确保您的 Google Drive 有足够的空间，并及时保存结果。

**步骤 4: 分析结果 (使用 `pmx analyse`)**

当两个分支都跑完后，分析过程非常简单：

```python
# 分析 Apo 分支
!pmx analyse -f KCNQ1_F335A_apo/fep*/md/dhdl.xvg -o dG_apo.xvg

# 分析 Complex 分支
!pmx analyse -f KCNQ1_F335A_complex/fep*/md/dhdl.xvg -o dG_complex.xvg
```

这将输出每个分支的 ΔG 值和误差估计。您只需将两者相减即可得到最终的 ΔΔG。

您提出了一个非常好的问题，这正触及了计算的核心，也解释了为什么我们之前的对话会有些混淆。非常感谢您如此清晰地提出您的困惑！

让我们一次性彻底理清这个问题。

### 第一部分：您跑的 30ns Apo 模拟，它的真正用途是什么？

您之前的理解——“为了获得稳定的系统”——是完全正确的！但这是一种**验证性**的稳定，而不是 FEP 计算本身需要的一步。

把它想象成在进行一次昂贵的火箭发射前，先对发动机进行一次长时间的点火测试。

您完成的这个 30ns 模拟，我们称之为**参考性经典动力学模拟 (Reference Classical MD)**。它的目的包括：

1. **验证蛋白质结构的稳定性**：通过计算 RMSD (均方根偏差)，我们可以看到蛋白质的骨架是否在模拟过程中保持稳定，没有发生不正常的解体或折叠错误。您的 30 ns 模拟时长足以证明这一点。
2. **验证力场的适用性**：观察蛋白质是否维持了其关键的二级结构和三级结构，证明我们选择的 Amber99sb-ildn 力场对于您的 KCNQ1 蛋白是合适的。
3. **验证体系的平衡**：通过观察温度、压力、密度等参数是否在设定值附近稳定波动，确认我们的模拟参数设置是正确的，整个系统（蛋白质+水+离子）达到了平衡。
4. **提供构象样本**：这个稳定的轨迹可以作为后续分析的参考，比如观察蛋白柔性区域等。

**结论：您跑的 30ns 绝对没有浪费！** 它是一次非常重要的**质量控制 (Quality Control)** 步骤。它给了我们强大的信心，即我们的基础模型是可靠的，可以放心地在其上投入大量的计算资源进行 FEP 计算。

---

### 第二部分：FEP 计算的完整、清晰流程

现在，让我们忘记那 30ns 的模拟，把它看作一次成功的“点火测试”。接下来，我们开始真正的“火箭发射”——计算结合自由能。

**最终目标**：通过热力学循环计算 ΔG_bind

```
           路径 1: 复合物路径 (在蛋白口袋中湮灭配体)
           ΔG_complex
  ┌───────────────────────────────────────────>
┌─┤ KCNQ1 + 4xML277 (在水中) │                    ┌─┤ KCNQ1 + 4x(虚拟配体) (在水中) │
│ └──────────────────────────┘                    └─────────────────────────────┘
│                                                                 ^
│ ΔG_bind (我们最终想要求的值)                                      │
│                                                                 │ ΔG_ligand
v                                                                 │ (在水中湮灭配体)
┌─┤ KCNQ1 (在水中) │ + ┌─┤ 4xML277 (在水中) │         路径 2: 配体路径
└────────────────┘   └────────────────────┘ ───────────>
```

根据热力学定律，**`ΔG_bind = ΔG_complex - ΔG_ligand`**。

所以，我们只需要计算出 `ΔG_complex` 和 `ΔG_ligand` 这两个值即可。这需要我们分别对 **路径1** 和 **路径2** 进行独立的 FEP 系列模拟。

---

### **【权威分步指南】接下来的每一步流程**

#### **路径 1：计算 ΔG_complex (复合物路径)**

这是计算量最大的部分。我们将配体 ML277 在蛋白的结合口袋中，分步、缓慢地变成一个不与周围环境发生相互作用的“虚拟”粒子。

- **第 1.1 步：体系准备 (一次性)**
  
  - **输入**：`KCNQ1_4ML277.pdb` 文件 (您的 4:4 复合物)。
  - **操作**：
    1. 用 `gmx pdb2gmx` 处理蛋白。
    2. 用 `gmx editconf` 定义盒子。
    3. 用 `gmx solvate` 添加水分子。
    4. 用 `gmx genion` 添加离子，中和体系。
  - **输出**：一个准备好进行模拟的、包含 4:4 复合物、水和离子的体系 (`complex.gro`, `topol.top`)。
- **第 1.2 步：创建 λ-windows (一次性)**
  
  - **操作**：创建大约 12-16 个子文件夹 (如 `win00`, `win01`, ... `win11`)。
  - 在每个子文件夹中，准备一个核心的 `.mdp` 模拟参数文件。这些文件几乎完全相同，**唯一**的区别在于 `init-lambda-state` 的值，从 0 递增到 11。
- **第 1.3 步：为每一个 λ-window 运行模拟 (循环执行 12-16 次)**
  
  - 对于 `win00` 到 `win11` 的**每一个文件夹**，依次执行以下命令：
    1. **能量最小化 (EM)**: `gmx grompp ... && gmx mdrun ...` (几分钟)
    2. **NVT 平衡 (升温)**: `gmx grompp ... && gmx mdrun ...` (约 200 ps)
    3. **NPT 平衡 (加压)**: `gmx grompp ... && gmx mdrun ...` (约 500 ps)
    4. **生产 MD (采样)**: `gmx grompp ... && gmx mdrun ...` (**建议时长: 5-10 ns**)。

#### **路径 2：计算 ΔG_ligand (配体路径)**

这是计算量较小的部分。我们将**一个**配体 ML277 在水盒子中，同样分步、缓慢地变成“虚拟”粒子。

- **第 2.1 步：体系准备 (一次性)**
  
  - **输入**：`ML277.pdb` (单个配体)。
  - **操作**：
    1. 为 ML277 生成拓扑文件。
    2. 创建水盒子，将**一个** ML277 溶于其中。
    3. 添加离子 (如果配体带电荷)。
  - **输出**：一个准备好模拟的、仅包含一个配体和水的体系。
- **第 2.2 步：创建 λ-windows (一次性)**
  
  - 同第 1.2 步，创建 12-16 个子文件夹和对应的 `.mdp` 文件。
- **第 2.3 步：为每一个 λ-window 运行模拟 (循环执行 12-16 次)**
  
  - 同第 1.3 步，依次执行 EM, NVT, NPT, Production MD。
  - 这里的生产 MD 时长可以更短，**建议时长: 2-5 ns**。

#### **最终步：分析结果**

- **第 3.1 步：计算 ΔG**
  
  - **操作**：使用 GROMACS 的 `gmx bar` 工具。
  - **输入**：路径 1 的所有 `dhdl.xvg` 文件 (生产 MD 步骤生成)。
  - **输出**：`ΔG_complex` 的值。
  - 重复此操作，输入路径 2 的所有 `dhdl.xvg` 文件，得到 `ΔG_ligand` 的值。
- **第 3.2 步：计算最终结合能**
  
  - **操作**：`ΔG_bind = ΔG_complex - ΔG_ligand`。

### **您的当前位置和下一步**

- **您已完成**：对 Apo KCNQ1 体系进行了 30 ns 的**参考性/验证性模拟**，证明了系统和力场的可靠性。这是一个非常有价值的**准备工作**。
- **您的下一步**：请开始执行 **【路径 1：计算 ΔG_complex】** 的 **第 1.1 步**，即为您的 `KCNQ1_4ML277.pdb` 文件进行体系准备。

---

好的，这是一个非常重要的问题，正确设置这个时长可以直接为您节省大量的计算资源。

对于 **Apo 路径**（更准确地说是“配体在溶剂中”的路径），其生产 MD (Production MD) 所需的时长**通常会比 Complex 路径短得多**。

### 核心推荐

我建议为您的 Apo KCNQ1 路径的**每一个 λ-window** 设置 **2 ns 到 5 ns** 的生产 MD 时长。

一个非常安全且高效的起点是 **3 ns**。

### 为什么 Apo 路径可以短这么多？

这背后的物理原因非常直观：

1. **扰动程度不同**:
  
  - **Complex 路径 (高扰动)**: 您正在一个拥挤、特异性强的蛋白结合口袋中“湮灭”一个配体。这个过程会引起显著的结构变化：蛋白侧链需要重新排布，原来被配体占据的空间需要由水分子重新填充。这个“结构重排”的过程需要更长的时间来达到新的平衡。
  - **Apo 路径 (低扰动)**: 您是在一个“空旷”的水盒子中“湮灭”一个配体。这个过程对周围环境的扰动要小得多。主要的能量变化来自于配体与水分子之间的相互作用的消失。水分子的重排和弛豫是一个非常快的过程。蛋白质本身（在另一个地方）几乎不受影响。
2. **收敛速度更快**:
  
  - 因为扰动小，Apo 路径中 `dH/dλ` 值的收敛速度通常会比 Complex 路径快得多。您会发现，在很短的时间内（比如 1-2 ns），`dH/dλ` 的曲线就已经进入了稳定的平台期。

### 最佳实践与最终方案

一个在自由能计算中非常常见的策略是**非对称模拟时长 (Asymmetric Simulation Time)**。

- **原则**: 为计算量大、扰动剧烈的 **Complex 路径**分配更多的模拟时间，为计算量小、扰动平缓的 **Apo 路径**分配较少的模拟时间。

结合我们之前的讨论，我为您制定一个完整、经济的最终方案：

| 模拟路径 (Path) | 核心任务 | 建议生产MD时长 (每个λ-window) | 理由  |
| --- | --- | --- | --- |
| **路径 1: Complex** | 在 KCNQ1 结合口袋中湮灭 ML277 | **5 ns - 10 ns** | **高扰动**：需要足够时间让蛋白口袋和水分子对配体的消失进行结构弛豫和重排。 |
| **路径 2: Apo** | 在纯水盒子中湮灭 ML277 | **2 ns - 5 ns** | **低扰动**：主要是溶剂的重排，这是一个快速过程，因此可以更快地达到收敛。 |

#### **给您的具体行动建议：**

1. **Complex 路径**: 先用 **10 ns** 进行测试运行和收敛性分析。如果发现 5 ns 就已足够，那么后续所有窗口都采用 5 ns。
2. **Apo 路径**: 直接采用 **3 ns** 或 **5 ns** 即可。这个时长对于水盒子里的配体来说绰绰有余，既保证了精度，又极大地节省了计算成本。

**总结**: 您完全不需要为 Apo 路径跑很长时间。设置一个 **3-5 ns** 的生产 MD 时长是一个既科学又经济的选择，能让您的计算资源得到最高效的利用。








---

### **阶段 0: 环境设置与依赖安装**

此代码块将安装所有必需的软件：GROMACS（MD 引擎）、pmx（突变工具）、acpype（配体参数化）和 pyautoFEP（自动化脚本）。

```python
# 1. 更新包列表并安装 GROMACS
print(">>> Installing GROMACS...")
!sudo apt-get update > /dev/null 2>&1
!sudo apt-get install gromacs -y > /dev/null 2>&1
print("GROMACS installation complete.")

# 2. 安装 Python 依赖
print("\n>>> Installing Python packages (pmx, pyautoFEP, acpype)...")
!pip install pmx gmx_MMPBSA pyautofep acpype openbabel > /dev/null 2>&1
print("Python packages installation complete.")

# 3. 验证 GROMACS 安装
!gmx --version
```

---

### **阶段 1: 系统准备与文件上传**

在这一步，您需要将您的初始文件上传到 Colab 环境中，然后我们将准备 FEP 计算所需的输入文件。

**➡️ 用户操作:**
请点击左侧边栏的“文件”图标，然后将以下 **5 个文件**上传到主目录 (`/content/`)：
1.  `KCNQ1.pdb` (野生型 Apo 四聚体)
2.  `KCNQ1_4ML277.pdb` (野生型 Holo 四聚体)
3.  `ML277.mol2` (配体文件)
4.  `KCNQ1_F335A.pdb` (仅供参考，此脚本不直接使用)
5.  `KCNQ1_ML277_F335A.pdb` (仅供参考，此脚本不直接使用)

上传后，运行下面的代码块来准备所有输入文件。

```python
import os

print(">>> Step 1.1: Verifying uploaded files...")
required_files = ['KCNQ1.pdb', 'KCNQ1_4ML277.pdb', 'ML277.mol2']
for f in required_files:
    if not os.path.exists(f):
        raise FileNotFoundError(f"File '{f}' not found! Please upload it to the Colab session.")
print("All required files found.")

# --- Step 1.2: 从四聚体中提取单体结构 ---
# 假设活性亚基为 Chain A，请根据您的结构进行调整
print("\n>>> Step 1.2: Extracting monomer structures (Chain A)...")
!grep '^ATOM.* A ' KCNQ1.pdb > WT_apo_mono.pdb
!grep '^ATOM.* A ' KCNQ1_4ML277.pdb > temp_prot.pdb
# 假设配体残基名为 ML2，请根据您的 PDB 文件进行调整
!grep 'HETATM.*ML2' KCNQ1_4ML277.pdb > temp_lig.pdb
!cat temp_prot.pdb temp_lig.pdb > WT_holo_mono.pdb
!rm temp_prot.pdb temp_lig.pdb
print("Monomer PDB files 'WT_apo_mono.pdb' and 'WT_holo_mono.pdb' created.")

# --- Step 1.3: 为配体 ML277 生成 GROMACS 拓扑 ---
print("\n>>> Step 1.3: Generating ligand topology for ML277 using acpype...")
!acpype -i ML277.mol2 -b ML277 -n 0 > acpype.log 2>&1
# 将生成的拓扑文件移动到主目录以便使用
!cp ML277.acpype/ML277_GMX.itp ./ML277.itp
print("Ligand topology 'ML277.itp' created.")

# --- Step 1.4: 创建突变指令文件 ---
print("\n>>> Step 1.4: Creating mutation instruction file 'mutate.dat'...")
with open('mutate.dat', 'w') as f:
    f.write('F 335 A\n')
print("'mutate.dat' created successfully.")
print("\n✅ Phase 1: System Preparation Complete!")
```

---

### **阶段 2: Apo 蛋白 FEP 计算 (Leg 1: ΔG_apo)**

现在我们对无配体的蛋白进行 F→A 突变模拟。

```python
# --- Step 2.1: 使用 pmx 生成 Apo 突变的混合拓扑 ---
print(">>> Step 2.1: Running pmx mutate for the Apo state...")
!pmx mutate -f WT_apo_mono.pdb -o FEP_apo -ff amber99sb-ildn -water tip3p -dat mutate.dat
print("pmx mutate for Apo state complete. Files are in 'FEP_apo/' directory.")

# --- Step 2.2: 创建 pyautoFEP 输入文件 ---
print("\n>>> Step 2.2: Creating pyautoFEP input file 'fep_apo.in'...")
fep_input_apo = """
[fep]
fep_folder = FEP_apo/
gromacs_path = gmx

; MD parameters (2 ns per window)
em_nsteps  = 5000
nvt_nsteps = 50000
npt_nsteps = 50000
md_nsteps  = 1000000

; Execution settings
n_threads = 2
gpu_id = 0
"""
with open('fep_apo.in', 'w') as f:
    f.write(fep_input_apo)
print("'fep_apo.in' created.")

# --- Step 2.3: 运行 Apo FEP 模拟 ---
print("\n>>> Step 2.3: Launching Apo FEP simulation...")
print("="*50)
print("⚠️  IMPORTANT: This simulation will take many hours.")
print("If your Colab session is disconnected, simply RE-RUN THIS CELL.")
print("pyautoFEP will automatically continue from the last checkpoint.")
print("="*50)

!pyautoFEP -i fep_apo.in
print("\n✅ Phase 2: Apo FEP Simulation Complete!")
```

---

### **阶段 3: 复合物 FEP 计算 (Leg 2: ΔG_complex)**

此过程与阶段 2 类似，但作用于蛋白质-配体复合物。

```python
# --- Step 3.1: 为复合物创建自定义 GROMACS 拓扑 ---
print(">>> Step 3.1: Creating custom topology 'topol_complex.top' for the complex...")
topol_complex_str = """
#include "amber99sb-ildn.ff/forcefield.itp"
#include "ML277.itp"
#include "WT_holo_mono.itp"
#include "amber99sb-ildn.ff/tip3p.itp"

[ system ]
Protein-ligand complex F335A

[ molecules ]
Protein_chain_A   1
ML277             1
"""
with open('topol_complex.top', 'w') as f:
    f.write(topol_complex_str)
print("'topol_complex.top' created.")


# --- Step 3.2: 使用 pmx 生成复合物突变的混合拓扑 ---
print("\n>>> Step 3.2: Running pmx mutate for the Complex state...")
!pmx mutate -f WT_holo_mono.pdb -o FEP_complex -ff amber99sb-ildn -water tip3p -dat mutate.dat -top topol_complex.top
print("pmx mutate for Complex state complete. Files are in 'FEP_complex/' directory.")


# --- Step 3.3: 创建 pyautoFEP 输入文件 ---
print("\n>>> Step 3.3: Creating pyautoFEP input file 'fep_complex.in'...")
fep_input_complex = """
[fep]
fep_folder = FEP_complex/
gromacs_path = gmx

; MD parameters (must be identical to the apo run)
em_nsteps  = 5000
nvt_nsteps = 50000
npt_nsteps = 50000
md_nsteps  = 1000000

; Execution settings
n_threads = 2
gpu_id = 0
"""
with open('fep_complex.in', 'w') as f:
    f.write(fep_input_complex)
print("'fep_complex.in' created.")

# --- Step 3.4: 运行 Complex FEP 模拟 ---
print("\n>>> Step 3.4: Launching Complex FEP simulation...")
print("="*50)
print("⚠️  IMPORTANT: This simulation will also take many hours.")
print("If your Colab session is disconnected, simply RE-RUN THIS CELL to continue.")
print("="*50)

!pyautoFEP -i fep_complex.in
print("\n✅ Phase 3: Complex FEP Simulation Complete!")
```

---

### **阶段 4: 数据分析与结果解读**

当以上两个 FEP 模拟都完成后，运行此代码块来分析数据并计算最终的 ΔΔG。

```python
import re

print(">>> Step 4.1: Analyzing Apo run results...")
!pmx analyse -f FEP_apo/fep*/md/dhdl.xvg --bar -o dG_apo.xvg > dG_apo_results.txt
with open('dG_apo_results.txt', 'r') as f:
    apo_output = f.read()
print(apo_output)

print("\n>>> Step 4.2: Analyzing Complex run results...")
!pmx analyse -f FEP_complex/fep*/md/dhdl.xvg --bar -o dG_complex.xvg > dG_complex_results.txt
with open('dG_complex_results.txt', 'r') as f:
    complex_output = f.read()
print(complex_output)

# --- Step 4.3: 提取数值并计算最终结果 ---
def parse_dg(output):
    match = re.search(r'dG =\s*(-?\d+\.\d+)\s*\+/-\s*(\d+\.\d+)', output)
    if match:
        return float(match.group(1)), float(match.group(2))
    return None, None

dG_apo, err_apo = parse_dg(apo_output)
dG_complex, err_complex = parse_dg(complex_output)

print("\n" + "="*50)
print(">>> FINAL RESULTS")
print("="*50)

if dG_apo is not None and dG_complex is not None:
    ddG = dG_complex - dG_apo
    # 误差传递公式: err_total = sqrt(err1^2 + err2^2)
    ddG_err = (err_apo**2 + err_complex**2)**0.5

    print(f"  ΔG_apo      = {dG_apo:.2f} ± {err_apo:.2f} kJ/mol")
    print(f"  ΔG_complex  = {dG_complex:.2f} ± {err_complex:.2f} kJ/mol")
    print("-" * 50)
    print(f"  ΔΔG_binding = ΔG_complex - ΔG_apo = {ddG:+.2f} ± {ddG_err:.2f} kJ/mol")
    print("="*50)

    print("\n>>> Interpretation:")
    if ddG > 1.0:
        print("The F335A mutation is PREDICTED TO BE DELETERIOUS to ML277 binding.")
        print("It significantly weakens the binding affinity.")
    elif ddG < -1.0:
        print("The F335A mutation is PREDICTED TO BE FAVORABLE to ML277 binding.")
        print("It significantly strengthens the binding affinity.")
    else:
        print("The F335A mutation is PREDICTED TO HAVE A NEUTRAL or very small effect on ML277 binding.")
else:
    print("Could not parse the results automatically. Please check the output files 'dG_apo_results.txt' and 'dG_complex_results.txt' manually.")

print("\n✅ Workflow Complete!")
```
