# GROMACS 分子动力学模拟 on Google Colab 标准操作流程 (v2 - 完整版)

本指南提供了一个在 Google Colab 环境中，利用 GPU 资源高效运行 GROMACS 分子动力学模拟的完整、标准化的工作流程。它包含了所有必要的脚本和参数文件，特别针对 Colab Pro/Pro+ 用户优化，可充分利用后台执行的优势。

---

## 目录

1. [**前提条件**](#前提条件)
2. [**第一步：环境设置与准备**](#第一步环境设置与准备)
3. [**第二步：安装 GROMACS**](#第二步安装-gromacs)
4. [**第三步：体系准备与预平衡**](#第三步体系准备与预平衡)
5. [**第四步：运行生产模拟 (核心)**](#第四步运行生产模拟-核心)
6. [**第五步：结果检查与基础分析**](#第五步结果检查与基础分析)

---

## 前提条件

在开始此流程前，请确保您已拥有：

- **初始结构文件**: 一个干净的蛋白质 `PDB` 文件 (例如: `protein.pdb`)。
- **Google Colab Pro 权限**: 强烈推荐，以便使用后台执行功能进行长时间模拟。

---

## 第一步：环境设置与准备

### 1.1 挂载 Google Drive

此步确保所有文件都能持久化保存，防止因 Colab 实例断开而丢失数据。

```python
from google.colab import drive
drive.mount('/content/drive')
```

### 1.2 创建并进入工作目录

在 Google Drive 中为项目创建专属文件夹，能让文件管理清晰有序。

```bash
# ‼️ 请根据您的实际情况修改为您的 Google Drive 路径
%cd /content/drive/MyDrive/FEP_F335A
!mkdir -p GMX_Simulation_v2
%cd GMX_Simulation_v2
```

> **笔记**: `%cd` 是 Colab 的 "magic command"，用于更改当前工作目录。后续所有命令都将在此目录下执行。

### 1.3 检查 GPU 分配

确认 Colab 已分配 GPU 资源。

```bash
!nvidia-smi
```

> **笔记**: 若未显示 GPU 信息，请前往 `代码执行程序` -> `更改运行时类型`，在 `硬件加速器` 中选择 `GPU`。

---

## 第二步：安装 GROMACS

使用 `Conda` 是在 Colab 中安装 GROMACS 最稳定、最推荐的方法。

```python
# 安装并初始化 Conda for Colab
!pip install -q condacolab
import condacolab
condacolab.install()
```

> **笔记**: 运行此单元格后，Colab 内核会自动重启，这是正常现象，请等待其完成。

```bash
# 在新的 Conda 环境中安装 GROMACS (指定版本以保证可复现性)
!mamba install -c conda-forge gromacs=2023.3
```

### 2.1 验证安装

检查 GROMACS 是否成功安装并能调用 GPU。

```bash
# 检查 gmx mdrun 是否能检测到 GPU (关键！)
!gmx mdrun -version
```

> **笔记**: 输出中应包含 `GPU support: CUDA` 字样，并列出检测到的 GPU 设备。

---

## 第三步：体系准备与预平衡

### 3.1 创建 MDP 参数文件

此步骤使用 `%%writefile` 命令自动生成所有需要的 `.mdp` 文件。**直接运行以下四个代码块即可。**

```bash
%%writefile minim.mdp
; 用于能量最小化的mdp文件
integrator  = steep         ; 使用最速下降法
emtol       = 1000.0        ; 能量最小化的收敛阈值 (kJ/mol/nm)
emstep      = 0.01          ; 初始步长
nsteps      = 50000         ; 最大步数
nstenergy   = 1000          ; 每1000步记录一次能量
nstlog      = 1000          ; 每1000步在log文件中记录一次
nstxout-compressed = 1000   ; 每1000步保存一次坐标
cutoff-scheme   = Verlet
coulombtype     = PME       ; 使用PME处理长程静电
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
constraints     = none      ; 能量最小化时不使用约束
```

```bash
%%writefile nvt.mdp
; 用于NVT平衡的mdp文件
title                   = NVT equilibration
define                  = -DPOSRES  ; 位置限制
integrator              = md        ; 蛙跳式积分器
nsteps                  = 50000     ; 模拟步数 (50000 * 2fs = 100ps)
dt                      = 0.002     ; 时间步长 (2 fs)
nstxout-compressed      = 5000      ; 每10ps保存一次坐标
nstlog                  = 5000      ; 每10ps记录一次log
nstenergy               = 5000      ; 每10ps记录一次能量
continuation            = no        ; 非延续模拟
constraint_algorithm    = lincs     ; 约束算法
constraints             = h-bonds   ; 约束与氢原子相连的键
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
tcoupl                  = V-rescale             ; 温度耦合方法
tc-grps                 = Protein Non-Protein   ; 分别控温的组
tau_t                   = 0.1   0.1             ; 时间常数 (ps)
ref_t                   = 300   300             ; 参考温度 (K)
pcoupl                  = no        ; NVT阶段不使用压力耦合
gen_vel                 = yes       ; 初始随机产生速度
gen_temp                = 300       ; 随机速度对应的温度
gen_seed                = -1        ; 随机数种子
```

```bash
%%writefile npt.mdp
; 用于NPT平衡的mdp文件
title                   = NPT equilibration
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000     ; 100 ps
dt                      = 0.002
nstxout-compressed      = 5000
nstlog                  = 5000
nstenergy               = 5000
continuation            = yes       ; 从NVT阶段延续
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1   0.1
ref_t                   = 300   300
pcoupl                  = Parrinello-Rahman     ; 压力耦合方法
pcoupltype              = isotropic             ; 各向同性耦合
tau_p                   = 2.0                   ; 时间常数 (ps)
ref_p                   = 1.0                   ; 参考压力 (bar)
compressibility         = 4.5e-5                ; 水的压缩率
refcoord_scaling        = com
```

```bash
%%writefile md.mdp
; 用于生产模拟的mdp文件
title                   = Production MD
integrator              = md
nsteps                  = 15000000  ; 15,000,000 * 2fs = 30,000 ps = 30 ns
dt                      = 0.002
nstxout-compressed      = 50000     ; 每100ps保存一次轨迹
nstlog                  = 50000     ; 每100ps记录一次log
nstenergy               = 50000     ; 每100ps记录一次能量
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1   0.1
ref_t                   = 300   300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
```

### 3.2 运行预处理流程

> **假设**: 您已将 `protein.pdb` 文件上传到当前工作目录。

```bash
# 1. 生成拓扑 (力场和水模型请根据体系选择)
# ‼️ 注意: -ff 和 -water 是关键参数，需根据您的分子和研究目的选择
!gmx pdb2gmx -f protein.pdb -o processed.gro -water tip3p -ff charmm27 -ignh

# 2. 定义盒子
!gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic

# 3. 溶剂化 (加水)
!gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# 4. 添加离子 - 生成tpr
!gmx grompp -f minim.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1

# 5. 添加离子 - 替换溶剂分子 (自动中和体系)
# ‼️ 注意: echo "SOL" 是为了自动选择要替换的组 (水分子)
!echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# 6. 能量最小化
!gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
!gmx mdrun -v -deffnm em

# 7. NVT 平衡
!gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
!gmx mdrun -v -deffnm nvt

# 8. NPT 平衡
!gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
!gmx mdrun -v -deffnm npt
```

---

## 第四步：运行生产模拟 (核心)

### 4.1 准备生产模拟 TPR 文件

使用 NPT 平衡的最终结构 `npt.gro` 和状态 `npt.cpt`。

```bash
!gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
```

### 4.2 **执行生产模拟 (可关闭浏览器)**

此命令经过优化，可充分利用 Colab 的 GPU 资源，并支持后台运行。

```bash
!gmx mdrun -v -deffnm md -nb gpu -pme gpu -nt 2
```

> #### **Colab Pro 用户须知**:
> 
> > **启动此命令后，您可以完全关闭浏览器或电脑。** 计算任务会在 Google 的云端服务器上持续进行。在预计的完成时间之后，再回到这个 Colab 笔记本，即可看到完整的日志输出和所有结果文件。

---

## 第五步：结果检查与基础分析

模拟结束后，回到您的 Colab 笔记本。

### 5.1 检查日志

检查 `md.log` 文件末尾，确保模拟正常结束。

```bash
!tail -n 30 md.log
```

> **笔记**: 您应看到 `Finished mdrun` 和最终的性能统计信息，代表模拟已成功完成。

### 5.2 基础轨迹分析 (示例)

处理周期性边界条件，并将蛋白质居中，然后计算 RMSD。

```bash
# 交互式选择 1 (Protein) 用于居中, 0 (System) 用于输出
!echo "1 0" | gmx trjconv -s md.tpr -f md.xtc -o md_final.xtc -pbc mol -center

# 交互式选择 4 (Backbone) 用于拟合, 4 (Backbone) 用于计算 RMSD
!echo "4 4" | gmx rms -s md.tpr -f md_final.xtc -o rmsd.xvg -tu ns
```

> **笔记**: 生成的 `rmsd.xvg` 文件可以下载到本地，使用 Grace, Xmgrace, 或 Python (Matplotlib) 等工具进行绘图。
