# 上传到 GitHub：纯网页操作版

本目录 `I:\work output\ChemWorkbench-github` 已经是**可以直接拖上去的仓库内容**：83 个文件、共约 850 KB，最大单文件 122 KB，远低于 GitHub 网页上传的 25 MB 限制。里面是 1.0.1 的源码和文档，**不含**打包好的 exe（源码包本来就没有）。

全程只需浏览器，不用装 Git。

**当前状态**：仓库 `https://github.com/paider-nk/Little-Chem-Bug` 已建好，Release `Little Chem 1.0.0` 已发布并挂上 115 MB 的 zip。仓库里目前只有 GitHub 自动生成的 README.md 和 LICENSE —— 第 2 步把那 83 个文件全部上传后，这两个同名文件会被直接覆盖成正式版本，不需要先删。

---

## 第 1 步：新建空仓库

1. 打开 <https://github.com>，登录。
2. 点右上角 **+** → **New repository**。
3. 填写：

| 字段 | 填什么 |
| --- | --- |
| Repository name | `Little-Chem-Bug`（已建好） |
| Description | 粘贴下面那段英文（也可留空，之后再改） |
| 可见性 | **Public** |
| Add a README file | **不勾** |
| Add .gitignore | **None** |
| Choose a license | **None** |

4. 点绿色 **Create repository**。

> 三个"Add …"都不要勾。这个目录里已经有 README.md、.gitignore 和 LICENSE，勾了会让远端多出一个初始提交，反而麻烦。

英文 Description（可选）：

```
Windows GUI for Gaussian and ORCA: build 3D inputs from SMILES, names, polymers, metal complexes or CIF files, batch-run jobs, then analyze results - log reports, HOMO/LUMO levels, ESP surfaces (Multiwfn + VMD), binding and solvation energies, density differences. Self-contained x64 build, no Python; offline input generation.
```

---

## 第 2 步：上传全部文件

1. 在新仓库页面点 **Add file ▾** → **Upload files**。
2. 打开文件资源管理器进入 `I:\work output\ChemWorkbench-github`。
3. **Ctrl+A 全选**，把这批文件和三个文件夹（`assets`、`examples`、`tests`）一起拖进浏览器的虚线框。

   注意：
   - 拖的是文件夹**里面的内容**，不要把 `ChemWorkbench-github` 这个外层文件夹本身拖进去（否则所有文件会多套一层目录）。
   - `.gitignore`、`.gitattributes` 是点开头的文件，在 Windows 里也照常可见；`Ctrl+A` 会包含它们，别漏掉。
   - 不要拖 `I:\化小区\ChemWorkbench-1.0.0-Windows-x64\` 那个解压目录——里面有 2345 个文件、单个 dll 超过 19 MB，既会被网页上传拦住，也绝不该进仓库。

4. 等进度条走完，页面下方确认列表里大约 83 个文件、没有任何错误提示。
5. 提交信息填 `ChemWorkbench 1.0.1 源码`，点绿色 **Commit changes**。

上传完成的仓库首页应该能看到 README.md 里渲染出的截图和中文说明。

---

## 第 3 步：占位符已填好（无需操作）

README.md 的 Releases 链接和 LICENSE 署名都已填成 `paider-nk`。仓库地址是 `https://github.com/paider-nk/Little-Chem-Bug`。以后想改成真名或换账号，在网页上点文件右上角的铅笔图标修改即可。

在仓库页面按 `t` 可以按文件名搜索；按 `.` 会打开网页版 VS Code，也能直接改文件。

---

## 第 4 步：发布可执行文件（Release）

仓库里只有源码，别人想用还得自己装 Python 打包。放一个 Release 附件就是"下载即用"。

✅ **这一步已经做完了**：你已发布 `Little Chem 1.0.0` 并挂上了 115 MB 的 zip，可以跳过。下面保留，是以后发 1.0.1 时用的（源码是 1.0.1，现有附件是 1.0.0）。

**方案 A（推荐）：先打 1.0.1 的包再发**

需要在装有 Python 3.13 的电脑上构建：

```powershell
cd "I:\化小区\ChemWorkbench-1.0.1-Source"
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt -r requirements-build.txt
.\.venv\Scripts\python.exe build_release.py
```

运行结束后会在 `build-release\` 下得到 `ChemWorkbench-1.0.1-Windows-x64.zip`，再按下面步骤传上去。

**方案 B：先把 1.0.0 传上去，当"历史版本"**

直接传现有的 `ChemWorkbench-1.0.0-Windows-x64.zip`，tag 写 `v1.0.0`。README 里的下载链接用的是 `releases/latest`，会指向最新那个 Release；之后补上 1.0.1 就会自动切过去。

### 创建 Release 的步骤

1. 仓库首页右侧 **Releases** → **Create a new release**（或直接访问 `https://github.com/paider-nk/Little-Chem-Bug/releases/new`）。
2. **Choose a tag** 填 `v1.0.1`（或 `v1.0.0`）→ 点 **Create new tag on publish**。
3. Title 填 `ChemWorkbench 1.0.1`。
4. 说明里写主要功能、支持的 Gaussian/ORCA 版本、必须完整解压并保留 `_internal` 目录。
5. **Attach binaries** 处拖入那个 zip，等它上传完（115 MB 需要几分钟）。
6. 点 **Publish release**。

单文件 100 MB 是"仓库文件"的限制；Release 附件单文件上限 2 GB，115 MB 完全没问题。

---

## 第 5 步（可选）：完善仓库首页

- 仓库首页右侧 **About** 齿轮 → 粘贴上面的英文 Description、Topics 填：`computational-chemistry`、`gaussian`、`orca`、`quantum-chemistry`、`rdkit`、`cheminformatics`、`molecular-modeling`、`esp-surface`、`multiwfn`、`vmd`、`windows`。
- 想换或加截图：上传图片到 `assets/`，然后在 README 里用 `![说明](assets/文件名.png)` 引用。

---

## 常见报错

| 现象 | 原因与处理 |
| --- | --- |
| `Yowza, that's a big file` | 单文件超过 25 MB（网页上传限制）。本目录没有这种文件，出现说明误拖了 exe/dll/zip。 |
| 提示 .gitignore 里的文件被忽略 | 网页上传**不读** .gitignore，可以无视；这条规则只对命令行 git 生效。 |
| 上传后中文文件名乱码 | 用 Chrome 或 Edge 上传，不要用手机浏览器。 |
| README 里的图片不显示 | 确认 `assets/screenshot-input.jpg` 已上传，且路径大小写一致。 |
| 想撤销一次网页提交 | 打开该 commit → 右上 `...` → **Revert**。 |

## 附：以后想改用命令行（可选）

网页上传适合首次发布；后续改动走 git 更省事。装好 <https://git-scm.com/download/win> 之后：

```powershell
cd "I:\work output\ChemWorkbench-github"
git init -b main
git add .
git commit -m "init"
git remote add origin https://github.com/paider-nk/Little-Chem-Bug.git
git push -u origin main
```

首次 push 会弹浏览器让你登录 GitHub（"Sign in with your browser"），不需要配 SSH key。
