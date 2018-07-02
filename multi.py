#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
import particlefilter as pf

D2R = math.pi / 180.0

plotinterval = 10
condition = 0

fig, ax = plt.subplots()
ax.set_xlabel("Y[m]")
ax.set_ylabel("X[m]")
ax.grid(True)
ax.axis('equal')

P = pf.ParticleFilter(1000)
P.scatter([0, -10, 0], [1., 2., 5.])
P2 = pf.ParticleFilter(1000)
P2.scatter([0, 10, 0], [1., 2., 5.])
landmark_change = -1
moving_vehicle = 1
moving_time = 100

while landmark_change < 10:
    landmark_change = landmark_change + 1
    if moving_vehicle == 0:
        moving_vehicle = 1
    else:
        moving_vehicle = 0

    if landmark_change >= 1:
        moving_time = 100
    for i in range(moving_time):
        if condition == 1:
            if moving_vehicle == 0:
                P.move([1., 0., 0.], [0.01, 0.01, 0.1], 1.)
            else:
                P2.move([1., 0., 0.], [0.01, 0.01, 0.1], 1.)
        else:
            if moving_vehicle == 0:
                P.move([1., 0., 0.], [0.01, 0.01, 0.1], 1.)
                if i % 10 == 0:
                    thetaLM = np.arctan2(P.y-P2.y, P.x-P2.x)/D2R
                    thetaML = np.arctan2(P2.y-P.y, P2.x-P.x)/D2R
                    dist = np.sqrt(np.power(P2.x - P.x, 2.0) + np.power(P2.y - P.y, 2.0))
                    P.observeLBL(dist, thetaML, thetaLM, 0.,[P2.x, P2.y, 0.], 0.5, 6.0, 3.0, 3)
                    P.resample()
            else:
                P2.move([1., 0., 0.], [0.01, 0.01, 0.1], 1.)
                if i % 10 == 0:
                    thetaML = np.arctan2(P.y-P2.y, P.x-P2.x)/D2R
                    thetaLM = np.arctan2(P2.y-P.y, P2.x-P.x)/D2R
                    dist = np.sqrt(np.power(P.x - P2.x, 2.0) + np.power(P.y - P2.y, 2.0))
                    P2.observeLBL(dist, thetaML, thetaLM, 0.,[P.x, P.y, 0.], 0.5, 6.0, 3.0, 3)
                    P2.resample()
        if i % plotinterval == 0:
                P.display(ax)
                P2.display(ax)
                plt.pause(.1)
plt.show()

