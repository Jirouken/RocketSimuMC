# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import g
from scipy.integrate import odeint
import pandas as pd
import pyproj
import simplekml
from quaternion import Quaternion
from mass import Mass
from force import Force
from wind import Wind
import angle_of_attack as aoa


class RocketSimulator(object):

    def __init__(self):
        self.T = 150
        self.dt = 0.01
        self.time = np.arange(0.0, self.T, self.dt)
        self.N = len(self.time)

    def initialize(self, design, launch_condition):
        """ 初期化 """
        self.name = design['name']
        self.m_af = design['m_af']
        self.I_af = design['I_af']
        self.CP = design['CP']
        self.CG_a = design['CG_a']
        self.d = design['d']
        self.area = np.pi * (self.d ** 2) / 4.0
        self.len_a = design['len_a']
        self.inertia_z0 = design['inertia_z0']
        self.inertia_zT = design['inertia_zT']
        self.engine = design['engine']
        self.me_total = design['me_total']
        self.me_prop = design['me_prop']
        self.len_e = design['len_e']
        self.d_e = design['d_e']

        self.p0 = np.array([0., 0., 0.])  # position(x, y, z)
        self.condition_name = launch_condition['name']
        self.theta0 = launch_condition['AngleOfFire']
        self.phi0 = launch_condition['azimuthal']
        self.launch_rod = launch_condition['launch_rod']
        self.v0 = np.array([0., 0., 0.])  # velocity(vx, vy, vz)
        self.ome0 = np.array([0., 0., 0.])
        self.density = launch_condition['density']
        self.wind_R = launch_condition['StandardWind']
        self.z_R = launch_condition['StandardHeight']
        self.beta = launch_condition['WindDirection']  # wind direction
        self.wind_direction = np.array([np.cos(self.beta), np.sin(self.beta), 0.0])
        self.qua_theta0 = np.array([np.cos(0.5 * self.theta0), np.sin(0.5 * self.theta0), 0., 0.])  # x軸theta[rad]回転, 射角
        self.qua_phi0 = np.array([np.cos(0.5 * self.phi0), 0., 0., np.sin(0.5 * self.phi0)])  # z軸phi[rad]回転, 方位角
        self.wind_direction = np.array([np.cos(self.beta), np.sin(self.beta), 0.0])

        self.engine_data = np.loadtxt(self.engine)

        self.force = Force(self.area, self.engine_data, self.T, self.density)
        self.thrust = self.force.thrust()

        self.mass = Mass(self.m_af, self.I_af, self.CG_a, self.len_a, self.inertia_z0, self.inertia_zT, self.me_total,
                         self.me_prop, self.len_e, self.d_e, self.force.burn_time, self.T)
        self.M = self.mass.mass()
        self.Me = self.mass.me_t()
        self.Me_dot = self.mass.me_dot()
        self.CG = self.mass.CG()
        self.CG_dot = self.mass.CG_dot()
        self.Ie = self.mass.iexg()
        self.Inertia = self.mass.inertia()
        self.Inertia_z = self.mass.inertia_z()
        self.Inertia_dot = self.mass.inertia_dot()
        self.Inertia_z_dot = self.mass.inertia_z_dot()

        self.wind = Wind(self.z_R, self.wind_R)

    def deriv_mc(self, pi, vi, quai, omei, t, nw, nb):
        """ 運動方程式 """
        qt = Quaternion()
        # 機軸座標系の推力方向ベクトル
        r_Ta = np.array([0., 0., 1.0])
        # 慣性座標系重力加速度
        gra = np.array([0., 0., -g])
        # 機軸座標系の空力中心位置
        r = np.array([0., 0., self.CG(t) - self.CP])
        # 慣性座標系の推力方向ベクトル
        r_T = qt.rotation(r_Ta, qt.coquat(quai))
        r_T /= np.linalg.norm(r_T)
        # 慣性テンソル
        I = np.diag([self.Inertia(t), self.Inertia(t), self.Inertia_z(t)])
        # 慣性テンソルの時間微分
        I_dot = np.diag([self.Inertia_dot(t), self.Inertia_dot(t), self.Inertia_z_dot(t)])
        # 慣性座標系対気速度
        beta = self.beta + nb
        v_air = (1 + nw) * self.wind.wind(pi[2]) * np.array([np.cos(beta), np.sin(beta), 0.0]) - vi
        # 迎角
        alpha = aoa.aoa(qt.rotation(v_air, quai))
        # ランチロッド垂直抗力
        N = 0
        # ランチロッド進行中
        if np.linalg.norm(pi) <= self.launch_rod and r_T[2] >= 0:
            Mg_ = self.M(t) * gra - np.dot(self.M(t) * gra, r_T) * r_T
            D_ = self.force.drag(alpha, v_air) - np.dot(self.force.drag(alpha, v_air), r_T) * r_T
            N = -Mg_ - D_
        # 慣性座標系加速度
        v_dot = gra + (self.thrust(t) * r_T + self.force.drag(alpha, v_air) + N) / self.M(t)
        # クォータニオンの導関数
        qua_dot = qt.qua_dot(omei, quai)
        # 機軸座標系角加速度
        ome_dot = np.linalg.solve(I, - np.cross(r, qt.rotation(self.force.drag(alpha, v_air), quai))
                                  - np.dot(I_dot, omei) - np.cross(omei, np.dot(I, omei)))
        # ランチロッド進行中
        if np.linalg.norm(pi) <= self.launch_rod:
            # ランチロッド進行中は姿勢が一定なので角加速度0とする
            ome_dot = np.array([0., 0., 0.])

        return vi, v_dot, qua_dot, ome_dot

    def simulate_mc(self, sd_w=0.1, sd_b=0.1):
        noise_w = np.random.randn(self.N) * sd_w
        noise_b = np.random.randn(self.N) * sd_b
        """ 数値計算 """
        qt = Quaternion()
        p = np.empty((self.N + 1, 3))
        v = np.empty((self.N + 1, 3))
        qua = np.empty((self.N + 1, 4))
        ome = np.empty((self.N + 1, 3))
        p[0] = self.p0
        v[0] = self.v0
        qua[0] = qt.product(self.qua_phi0, self.qua_theta0)
        ome[0] = self.ome0
        count = 0

        for (i, t) in enumerate(self.time):
            # Euler method
            p_dot, v_dot, qua_dot, ome_dot = self.deriv_mc(p[i], v[i], qua[i], ome[i], t, noise_w[i], noise_b[i])

            # if np.isnan(qua_dot).any() or np.isinf(qua_dot).any() or np.isnan(ome_dot).any() or np.isinf(ome_dot).any():
            #     count = i
            #     break

            p[i + 1] = p[i] + p_dot * self.dt
            v[i + 1] = v[i] + v_dot * self.dt
            qua[i + 1] = qua[i] + qua_dot * self.dt
            ome[i + 1] = ome[i] + ome_dot * self.dt

            # vz<0かつz<0のとき計算を中断
            if v[i + 1][2] < 0 and p[i + 1][2] < 0:
                count = i + 1
                break

            qua[i + 1] /= np.linalg.norm(qua[i + 1])

            if t <= self.force.burn_time:
                p[i + 1][2] = max(0., p[i + 1][2])

            # vz<0かつz<0のとき計算を中断
            if v[i + 1][2] < 0 and p[i + 1][2] < 0:
                count = i + 1
                break

        self.p = p[:count + 1]
        self.v = v[:count + 1]
        self.qua = qua[:count + 1]
        self.ome = ome[:count + 1]

    def falling_range(self, n=1000, sd_w=0.1, sd_b=0.1):
        falling_area = []
        max_alttitude = []
        self.paths = []
        for i in range(n):
            self.simulate_mc(sd_w, sd_b)
            falling_area.append(self.p[-1])
            max_alttitude.append(self.p[:, 2].max())
            self.paths.append(self.p)
            if (i + 1) % 10 == 0:
                print(str((i+1)/n * 100) + '%')
        self.paths = np.array(self.paths)
        return np.array(falling_area), np.array(max_alttitude)

    def output_kml(self, place):
        # 原点からの距離
        def dis2d(x, y):
            return pow(pow(x, 2) + pow(y, 2), 0.5)

        # y軸とベクトル(x, y)のなす角
        def ang2d(x, y):
            # y軸上
            if x == 0:
                if y >= 0:
                    return 0.0
                else:
                    return 180.0
            # x軸上
            if y == 0:
                if x >= 0:
                    return 90.0
                else:
                    return 270.0
            # 第1象限
            if x > 0 and y > 0:
                return np.arctan(x / y) * 180 / np.pi
            # 第2象限
            if x > 0 and y < 0:
                return 180.0 + np.arctan(x / y) * 180 / np.pi
            # 第3象限
            if x < 0 and y < 0:
                return 180.0 + np.arctan(x / y) * 180 / np.pi
            # 第4象限
            if x < 0 and y > 0:
                return 360.0 + np.arctan(x / y) * 180 / np.pi

        distance = [dis2d(self.p[i, 0], self.p[i, 1]) for i in range(len(self.p))]

        angle = [ang2d(self.p[i, 0], self.p[i, 1]) for i in range(len(self.p))]

        coordinate0 = place[1]
        latitude = coordinate0[0]
        longitude = coordinate0[1]
        geod = pyproj.Geod(ellps='WGS84')
        newLong = []
        newLat = []
        invAngle = []
        for i in range(len(self.p)):
            nlon, nlat, nang = geod.fwd(longitude, latitude, angle[i], distance[i])
            newLong.append(nlon)
            newLat.append(nlat)
            invAngle.append(nang)

        kml = simplekml.Kml(open=1)
        cood = []
        for i in range(len(self.p)):
            cood.append((newLong[i], newLat[i], self.p[i, 2]))
        ls = kml.newlinestring(name=self.name+"'s Path")
        ls.coords = cood
        ls.style.linestyle.width = 3
        ls.style.linestyle.color = simplekml.Color.blue
        ls.altitudemode = simplekml.AltitudeMode.relativetoground
        ls.extrude = 0
        kml.save(self.name + "_" + place[0] + '_' + self.condition_name + "_path.kml")
        print("Kml file for Google Earth was created.")


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from matplotlib import mlab
    from sklearn.mixture import GaussianMixture

    # 設計データ読み込み
    design = pd.read_csv('design.csv')
    # 打ち上げ条件読み込み
    condition = pd.read_csv('condition.csv')

    # モンテカルロシミュレーション
    rs = RocketSimulator()
    rs.initialize(design.loc[0], condition.loc[2])
    fr, ma = rs.falling_range(n=1000)

    # ガウシアンフィッティング
    n_components = 1
    gm = GaussianMixture(n_components, covariance_type='full')
    gm.fit(fr)

    xx = np.linspace(-7.0, 7.0, 200)
    yy = np.linspace(-5.0, 5.0, 200)
    X, Y = np.meshgrid(xx, yy)

    Z = mlab.bivariate_normal(X, Y, np.sqrt(gm.covariances_[0][0][0]), np.sqrt(gm.covariances_[0][1][1]),
                              gm.means_[0][0], gm.means_[0][1], gm.covariances_[0][0][1])

    plt.figure(figsize=(8, 8))
    plt.scatter(0, 0, label='launch complex', marker='*', c='r', s=100)
    plt.scatter(fr[:, 0], fr[:, 1], alpha=0.5, s=10)
    plt.xlim(-20, 20)
    plt.ylim(-90, 2)
    CS = plt.contour(X, Y, Z11)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('x', fontsize=20)
    plt.ylabel('y', fontsize=20)
