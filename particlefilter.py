# -*- coding: utf-8 -*-

import time
import math
import numpy as np
import matplotlib.pyplot as plt

D2R = math.pi / 180.0  # degree to radian

class ParticleFilter:
    def __init__(self, num):
        self.n = int(num)  # number of particles
        self.p = np.zeros((3, self.n))  # particles: horizontal position(x; y) [m], heading [deg]
        self.w = np.ones((1, self.n))  # weight
        self.x, self.y, self.a = 0., 0., 0.  # mean
        self.ex, self.ey, self.ea = 0., 0., 0.  # standard deviation

    def scatter(self, mean, sigma):  #パーティクルを散らす
        for i in range(3):
            self.p[i, :] = self.gen_rans(mean[i], sigma[i])
        self.w = np.ones((1, self.n))  # 重みは1でリセット
        self.gen_est()  #統計値の更新

    def gen_est(self):  # generate statistics values
        self.x, self.y, self.a = self.p.mean(1)
        self.ex, self.ey, self.ea = self.p.std(1)

    def display(self, axis):
        # plt.plot(P.w[0, :], '.')
        axis.plot(self.p[1,:], self.p[0,:], '.')
        #axis.show()

    def gen_rans(self, x, ex):  # generate random numbers with the length of self.n
        if ex <= 0:
            rans = np.zeros(self.n) + x
        else:
            rans = np.random.normal(x, ex, self.n)
        return rans

    def move(self, vin, ev, dt):
        # PREDICTION PHASE: move particles based on the motion measurements
        # v: velocity  ev: std. of velocity (u [m/s], v [m/s], r [deg/s])
        # dt: period [sec]
        v = np.zeros((3, self.n))
        for i in range(3):
            v[i, :] = self.gen_rans(vin[i], ev[i]) * dt
        dp = np.zeros((3, self.n))
        cosh = np.cos(self.p[2,:] * D2R)
        sinh = np.sin(self.p[2,:] * D2R)        
        for i in range(self.n):
            R = np.array([[cosh[i], -sinh[i]], [sinh[i], cosh[i]]])
            dp[0:2, i] = R.dot(v[0:2, i])
        dp[2, :] = v[2, :]
        self.p = self.p + dp
        self.gen_est()

    def resample(self):
        n_one = np.count_nonzero(self.w == 1)  # number of the samples with default weight
        if n_one < self.n:
            wp = np.concatenate((self.w, self.p))
            wphigh = wp[:, np.nonzero(self.w > 1.0)[1]]  # weight > 1 -> sorted
            whigh = wphigh[0,:].sum()
            n_pick = int(round(whigh / (whigh + n_one) * self.n))  # number of samples to be picked from the highly weighted group
            p2 = np.zeros((3, n_pick))
            i = 0
            k = 0
            wsum = wphigh[0, k]
            wcur = 0
            winc = whigh / float(n_pick)
            while i < n_pick:
                while wsum <= wcur:
                    k += 1
                    wsum += wphigh[0, k]
                p2[:,i] = wphigh[1:4, k]  # wphigh might be better to be sorted before resampling
                i += 1
                wcur += winc
            if n_pick < self.n:
                wpone = wp[:, np.nonzero(self.w == 1.0)[1]].T
                np.random.shuffle(wpone)
                wpone = wpone.T
                self.p = np.concatenate((p2, wpone[1:4, 0:(self.n - n_pick)]), 1)
            else:
                self.p = p2
            self.w = np.ones((1, self.n))
            # print wphigh.shape, n_pick
        self.gen_est()

    def observeLBL(self,range,thetaML, thetaLM, auv_z, pos, er, ethetaML, ethetaLM, k):
        dx = self.p[0, :] - pos[0]
        dy = self.p[1, :] - pos[1]
        dz = np.zeros((1, self.n)) + (auv_z - pos[2])
        dxyz = np.array([dx, dy, dz])
        r = np.linalg.norm(dxyz, axis=0)
        dr = r - range
        wtmp = np.exp(-dr**2/(2.0 * er**2)) * np.exp(k**2/2.0)

        dthetaML = np.arctan2(pos[1] - self.p[1,:], pos[0]-self.p[0,:]) / D2R - thetaML
        i=0
        while i < self.n:
            if dthetaML[i] > 180.0:
                dthetaML[i] = dthetaML[i] - 360.0
            if dthetaML[i] < -180.0:
                dthetaML[i] = dthetaML[i] +360.0
            i=i+1
        wtmp2 = np.exp(-dthetaML**2/(2.0 * ethetaML**2)) * np.exp(k**2/2.0)

        dthetaLM = np.arctan2(self.p[1,:] - pos[1], self.p[0,:] - pos[0]) /D2R - thetaLM
        i = 0
        while i < self.n:
            if dthetaLM[i] > 180.0:
                dthetaLM[i] = dthetaLM[i] - 360.0
            if dthetaLM[i] < -180.0:
                dthetaLM[i] = dthetaLM[i] + 360.0
            i = i+1
        wtmp3 = np.exp(-dthetaLM**2 / (2.0 * ethetaLM**2)) * np.exp(k** 2 / 2.0)
        wtmp = np.array([wtmp*wtmp2*wtmp3, np.ones((1, self.n))]).max(0)
        self.w = self.w * wtmp