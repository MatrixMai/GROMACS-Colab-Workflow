# **KCNQ1 F335A 突变对 ML277 结合能影响的 FEP 计算工作流**

### **目标**
使用自由能微扰（FEP）方法，计算 KCNQ1 蛋白的 F335A 点突变对配体 ML277 结合自由能的影响（ΔΔG）。

### **核心策略：热力学循环**
我们将通过计算两个独立的、计算上可行的“炼金术”突变路径的自由能（ΔG）来达成目标：
1.  **ΔG_apo**: 在无配体的 Apo 蛋白中，将 F335 突变为 A335。
2.  **ΔG_complex**: 在与 ML277 结合的 Holo 复合物中，进行同样的 F335A 突变。

最终的结合自由能变化为： **\[ ΔΔG_{binding} = ΔG_{complex} - ΔG_{apo} \]**

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
