# Mini Classical Compensator（迷你经典校正器）

**从零开始、零依赖的 C 语言实现**，涵盖经典校正器设计与 PID 控制理论。每个模块对应 MIT（及其他顶尖大学）的一门或多门课程，将教科书中的传递函数和频域设计方法转化为可运行的 C 代码，实现理论与实践的桥接。

## 模块状态：COMPLETE ✅

| 子模块 | include/+src/ 行数 | L1-L6 | L7 | L8 | L9 | 状态 |
|--------|--------------------|-------|----|----|----|------|
| mini-cascade-control | 3,194 | 已完成 | 已完成 | 部分 | 部分 | ✅ |
| mini-feedforward-control | 4,061 | 已完成 | 已完成 | 部分 | 部分 | ✅ |
| mini-lag-compensator | 5,362 | 已完成 | 已完成 | 部分 | 部分 | ✅ |
| mini-lead-compensator | 3,404 | 已完成 | 已完成 | 部分 | 部分 | ✅ |
| mini-lead-lag-design | 3,009 | 已完成 | 已完成 | 部分 | 部分 | ✅ |
| mini-pid-theory | 3,357 | 已完成 | 已完成 | 部分 | 部分 | ✅ |
| mini-pid-tuning-ziegler | 6,771 | 已完成 | 已完成 | 部分 | 部分 | ✅ |
| mini-ratio-control | 3,607 | 已完成 | 已完成 | 部分 | 部分 | ✅ |

全部 8 个子模块均超过 3,000 行代码阈值，`make test` 零失败通过。

## 子模块总览

| 子模块 | 主题 | 参考课程 |
|--------|------|----------|
| [mini-cascade-control](mini-cascade-control/) | 串级控制架构、嵌套回路、带宽分离、顺序闭环 | MIT 6.302, Stanford ENGR105, Berkeley ME232 |
| [mini-feedforward-control](mini-feedforward-control/) | 前馈补偿、二自由度架构、扰动抑制、输入整形 | MIT 6.302, Stanford ENGR105, Cambridge 3F2 |
| [mini-lag-compensator](mini-lag-compensator/) | 滞后校正器、稳态误差消除、直流增益提升、相位滞后分析 | MIT 6.302, Stanford ENGR105, Berkeley ME132 |
| [mini-lead-compensator](mini-lead-compensator/) | 超前校正器、相位裕度改善、瞬态响应、Bode 设计 | MIT 6.302, Stanford ENGR105, Caltech CDS 110 |
| [mini-lead-lag-design](mini-lead-lag-design/) | 超前-滞后综合、根轨迹设计、回路整形、频域方法、灵敏度 | MIT 6.302, Georgia Tech ECE 6550, ETH 151-0591 |
| [mini-pid-theory](mini-pid-theory/) | PID 形式、抗积分饱和、二自由度 PID、微分滤波、无扰切换 | MIT 6.302, Stanford ENGR105, Cambridge 3F2 |
| [mini-pid-tuning-ziegler](mini-pid-tuning-ziegler/) | Ziegler-Nichols 整定、FOPDT 辨识、继电自整定、Cohen-Coon、IMC | MIT 6.302, Chalmers (Åström & Hägglund) |
| [mini-ratio-control](mini-ratio-control/) | 比值控制、交叉限幅、混合控制、主从站、质量平衡完整性 | Purdue ECE 602, 工业过程控制 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块将经典控制传递函数与设计流程转化为可运行的算法
- **实用演示程序** — 校正器设计工具、PID 自整定器、串级整定、比值站仿真等

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-cascade-control
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-classical-compensator/
├── mini-cascade-control/        # 串级控制（嵌套主/副回路）
├── mini-feedforward-control/    # 前馈与二自由度控制架构
├── mini-lag-compensator/        # 相位滞后校正器（稳态精度）
├── mini-lead-compensator/       # 相位超前校正器（瞬态改善）
├── mini-lead-lag-design/        # 超前-滞后综合设计与回路整形
├── mini-pid-theory/             # PID 控制理论、形式与高级特性
├── mini-pid-tuning-ziegler/     # Ziegler-Nichols 及经典 PID 整定方法
└── mini-ratio-control/          # 比值控制（交叉限幅与混合）
```

## 许可证

MIT
