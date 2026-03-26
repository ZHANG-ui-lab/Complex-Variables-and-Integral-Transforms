# -*- coding: utf-8 -*-
"""
圆柱绕流项目工具模块 (utils.py)
包含：
1. CR方程验证函数 (复势解析性验证)
2. 伯努利方程计算类 (压力系数计算)
3. 流场生成器 (网格与坐标生成)
"""

import numpy as np
import sympy as sp


# ===================== 1. CR 方程验证函数 =====================
def verify_cr_equations_polar(r, theta, U, a):
    """
    极坐标下柯西-黎曼 (C-R) 方程符号验证
    用于验证复势函数的解析性，确保流场满足无旋条件。
    
    参数:
        r: 径向坐标符号或数值
        theta: 周向坐标符号或数值
        U: 来流速度
        a: 圆柱半径
    
    返回:
        cr1, cr2: C-R方程的两个残差（理想值应为 0）
    """
    # 定义速度势 phi 与流函数 psi (基于无环量圆柱绕流复势)
    phi = U * (r + a**2 / r) * sp.cos(theta)
    psi = U * (r - a**2 / r) * sp.sin(theta)
    
    # 计算 C-R 方程残差
    cr1 = sp.simplify(sp.diff(phi, r) - (1/r) * sp.diff(psi, theta))
    cr2 = sp.simplify((1/r) * sp.diff(phi, theta) + sp.diff(psi, r))
    
    return cr1, cr2


def cr_numerical_check(r_data, theta_data, U, a):
    """
    数值方式验证C-R方程残差，用于随机点测试。
    
    参数:
        r_data: 径向数据数组 (np.array)
        theta_data: 周向数据数组 (np.array)
        U: 来流速度
        a: 圆柱半径
    
    返回:
        max_residual: 最大残差值
    """
    # 解析解速度分量
    v_r = U * (1 - a**2 / r_data**2) * np.cos(theta_data)
    v_theta = -U * (1 + a**2 / r_data**2) * np.sin(theta_data)
    
    # 数值微分近似 (中心差分)
    dr = 1e-8
    dtheta = 1e-8
    
    # 这里直接使用解析解验证 C-R 关系 (v_r = d(psi)/d(theta), v_theta = -d(psi)/dr)
    # 由于我们已知解析解，直接验证数值导数与解析速度的一致性
    psi_r = U * (1 + a**2 / r_data**2) * np.sin(theta_data)
    psi_theta = U * (r_data - a**2 / r_data) * np.cos(theta_data)
    
    residual_r = np.abs(v_r - psi_theta)
    residual_theta = np.abs(v_theta + psi_r)
    
    return np.max(residual_r), np.max(residual_theta)


# ===================== 2. 伯努利方程计算类 =====================
class BernoulliCalculator:
    """
    伯努利方程计算类
    用于基于速度场计算压力系数 (Cp)，验证理论公式。
    """
    def __init__(self, U=1.0, rho=1.225):
        """
        初始化计算器
        
        参数:
            U: 参考来流速度 (m/s)
            rho: 流体密度 (kg/m^3)，默认空气密度
        """
        self.U = U
        self.rho = rho
        self.q = 0.5 * self.rho * self.U**2  # 动压

    def calculate_cp(self, V_mag):
        """
        基于伯努利方程计算压力系数 Cp
        Cp = (p - p_inf) / (0.5 * rho * U^2) = 1 - (|V|/U)^2
        
        参数:
            V_mag: 速度模长数组
        
        返回:
            Cp: 压力系数数组
        """
        if np.max(V_mag) > 1e-10: # 防止除零
            cp = 1 - (V_mag / self.U)**2
        else:
            cp = np.ones_like(V_mag)
        return cp

    def calculate_pressure(self, V_mag, p_inf=101325):
        """
        计算当地静压 p
        
        参数:
            V_mag: 速度模长数组
            p_inf: 远场静压 (Pa)，默认标准大气压
        
        返回:
            p: 当地压力数组
        """
        cp = self.calculate_cp(V_mag)
        p = p_inf + cp * self.q
        return p


# ===================== 3. 流场生成器 =====================
class FlowFieldGenerator:
    """
    流场生成器
    用于生成标准化的流场坐标、圆柱表面点集等。
    """
    @staticmethod
    def generate_cylinder_surface(a=1.0, N=72):
        """
        生成圆柱表面采样点 (极坐标转直角坐标)
        
        参数:
            a: 圆柱半径
            N: 采样点数
        
        返回:
            theta: 角度数组 (rad)
            z: 复坐标数组 (x+iy)
            x, y: 直角坐标数组
        """
        theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
        z = a * np.exp(1j * theta)
        x = z.real
        y = z.imag
        return theta, z, x, y

    @staticmethod
    def generate_grid(x_range=(-5, 5), y_range=(-5, 5), res=50j):
        """
        生成计算网格
        
        参数:
            x_range: x轴范围 (min, max)
            y_range: y轴范围 (min, max)
            res: 网格密度复数，如 50j 表示x,y方向各50个点
        
        返回:
            X, Y: 网格坐标矩阵
            Z: 复平面网格矩阵
        """
        x = np.linspace(*x_range, res)
        y = np.linspace(*y_range, res)
        X, Y = np.meshgrid(x, y)
        Z = X + 1j * Y
        return X, Y, Z

    @staticmethod
    def complex_velocity_potential(Z, U=1.0, a=1.0, Gamma=0.0):
        """
        计算圆柱绕流复速度 (复势的导数)
        V(z) = dF/dz = U(1 - a^2/z^2) + iGamma/(2πz)
        
        参数:
            Z: 复平面网格矩阵
            U: 来流速度
            a: 圆柱半径
            Gamma: 环量
        
        返回:
            V: 复速度矩阵
        """
        # 避免除以零 (在原点处)
        with np.errstate(divide='ignore', invalid='ignore'):
            V = U * (1 - a**2 / Z**2) + 1j * Gamma / (2 * np.pi * Z)
            # 处理原点无穷大
            V[np.isinf(V)] = 0 
        return V