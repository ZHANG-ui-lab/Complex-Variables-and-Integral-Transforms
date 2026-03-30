# Complex-Variables-and-Integral-Transforms
复变函数与积分变换课程项目 —— 圆柱绕流复势分析与可视化

## 项目简介
本项目基于**复变函数与势流理论**，实现理想不可压缩流体圆柱绕流的数学建模、数值计算与可视化分析。通过复势叠加原理构造流场，完成流场可视化、气动参数计算与误差验证，为后续多学科优化与环量修正提供基础支撑。

---

## 项目结构
```
Complex-Variables-and-Integral-Transforms/
├── 1_复势推导.ipynb          # 任务1：复势构造与解析性验证（C-R方程验证）
├── 2_流场可视化.ipynb         # 任务2.1：流场流线、等势线与速度场可视化
├── 3_压力分析.ipynb           # 任务2.2：圆柱表面压力系数计算与误差分析
├── utils.py                   # 工具模块（CR验证、伯努利计算、流场生成）
├── requirements.txt           # Python依赖库清单
└── README.md                  # 项目说明文档
```

---

## 核心功能
### 1. 复势模型与解析性验证
- 构造均匀来流 + 偶极子叠加的圆柱绕流复势
- 极坐标下**柯西-黎曼（C-R）方程**符号验证与数值校验
- 验证复势函数在流场区域的解析性，保证流场无旋、不可压缩

### 2. 流场可视化
- 生成流线、等势线分布，展示圆柱绕流的典型流态
- 绘制速度幅值云图，标注驻点位置
- 实现交互式可视化，直观展示流场特征

### 3. 气动参数计算与分析
- 基于伯努利方程计算圆柱表面**压力系数 \(C_p\)**
- 对比数值解与理论解 \(C_p = 1 - 4\sin^2\theta\)，验证计算精度（RMSE < 0.01）
- 生成极坐标/直角坐标系压力分布曲线，标注关键特征点
- 导出标准化 CSV 数据，支持后续分析与报告生成

---

## 快速开始
### 1. 环境配置
安装依赖库：
```bash
pip install -r requirements.txt
```
或手动安装：
```bash
pip install numpy sympy matplotlib pandas ipywidgets
```

### 2. 运行流程
1.  打开 `1_复势推导.ipynb`，完成复势构造与解析性验证
2.  运行 `2_流场可视化.ipynb`，生成流场可视化结果
3.  执行 `3_压力分析.ipynb`，完成压力系数计算与误差校验

### 3. 工具模块使用
```python
from utils import verify_cr_equations_polar, BernoulliCalculator, FlowFieldGenerator

# 生成圆柱表面采样点
generator = FlowFieldGenerator()
theta, z, x, y = generator.generate_cylinder_surface(N=72)

# 计算压力系数
calc = BernoulliCalculator(U=1.0)
V_mag = np.abs(generator.complex_velocity_potential(z))
cp = calc.calculate_cp(V_mag)
```

---

## 关键结果
- **解析性验证**：C-R方程残差接近 0，复势函数满足解析性要求
- **压力系数精度**：数值解与理论解 RMSE < 1e-15，远优于项目阈值 0.01
- **流场特征**：清晰展示驻点（\(C_p=1\)）、最小压力点（\(C_p=-3\)）等典型气动特征

---

## 技术栈
- **数学建模**：复变函数、势流理论、柯西-黎曼方程
- **数值计算**：NumPy、SymPy（符号计算）
- **可视化**：Matplotlib、ipywidgets（交互式）
- **开发环境**：Jupyter Notebook / VS Code

---

## 项目成员
- 张明悦
- 任孙亚舟
- 邓梓龙
- 张子琪
- 丁辰昊
- 李瑞烽

---

## 许可证
本项目仅用于课程学习与学术交流，未经允许不得用于商业用途。