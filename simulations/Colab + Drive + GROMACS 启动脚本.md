### **Colab + Drive + GROMACS 启动脚本**

这个脚本会自动完成所有必要的检查和设置：

1.  挂载您的Google Drive。
2.  检查在**当前会话**中，GROMACS是否已经被安装。
3.  如果没安装，它会自动找到您云端硬盘里的编译目录，并快速执行安装命令。
4.  自动将工作目录切换到您在Google Drive上的项目文件夹。
5.  最后验证安装是否成功。

```python
import os
from google.colab import drive

# 1. 挂载 Google Drive
print("▶️ 正在挂载 Google Drive...")
drive.mount('/content/drive')

# 2. 定义关键路径
# 你在Google Drive上的项目文件夹路径
gdrive_project_path = "/content/drive/MyDrive/FEP_F335A" 
# 你在Google Drive上GROMACS编译好的build目录路径
gromacs_build_path = os.path.join(gdrive_project_path, "gromacs-2023.3/build")
# GROMACS安装后，gmx命令的预期存放路径
gmx_executable_path = "/usr/local/gromacs/bin/gmx"

# 3. 检查GROMACS是否已在当前会话中安装
if os.path.exists(gmx_executable_path):
    print(f"✅ GROMACS 已安装在 {gmx_executable_path}，无需重复操作。")
else:
    print("▶️ 在当前会话中未找到 GROMACS，正在从您的云端硬盘恢复安装...")
    if not os.path.exists(gromacs_build_path):
        print(f"❌ 错误：在 {gromacs_build_path} 未找到 build 目录！")
        print("   请确保您已成功执行过一次完整的编译，并且文件已保存在云端硬盘。")
    else:
        # 进入build目录，并执行 'make install'
        # 这一步非常快，因为它只是复制已编译好的文件
        %cd {gromacs_build_path}
        !make install > /dev/null 2>&1
        print("✅ 'make install' 已完成。")

# 4. 进入你的项目工作目录，方便后续操作
print(f"\n▶️ 正在将当前目录切换至您的项目文件夹: {gdrive_project_path}")
%cd {gdrive_project_path}

# 5. 最终验证
print("\n" + "="*50)
print("最终验证:")
print(f"当前工作目录: ")
!pwd
print("\nGROMACS 版本信息:")
# 使用绝对路径确保执行的是我们刚刚安装的版本
!{gmx_executable_path} --version
print("="*50 + "\n")
print("✅ 环境设置完毕，您可以开始运行模拟了。")
```

### **核心要点总结**

| 项目 | 存储位置 | 持久性 |
| :--- | :--- | :--- |
| **您的数据和源代码** | Google Drive | **永久** |
| **GROMACS软件程序** | Colab 虚拟机 (`/usr/local/`) | **临时** |

通过每次运行这个启动脚本，您就实现了“即插即用”的体验，无需再为环境问题烦恼。
