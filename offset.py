#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def offset(x1, y1 ,x2, y2, xs, ys, theta):
    wx = x2 - x1
    wy = y2 - y1 
    x0 = np.cos(theta) 
    y0 = np.sin(theta) # 角からxの単位ベクトルを求む


    w = np.array([wx, wy]) #w1w2のなすベクトル
    x = np.array([x0, y0]) #auvの単位ベクトル

    wx = np.dot(w, x)
    wx_cross = np.cross(x,w) # 角の正負判定用の外積
    xl = np.linalg.norm(x)
    wl = np.linalg.norm(w)

    cos_fi = wx / (xl * wl) #内積よりオフセット角fiをもとめる
    fi = np.arccos(cos_fi)
    fi_deg = fi * 180 / 3.1415

    x_vector = np.array([x2 - xs, y2 - ys]) # AUVとW2を結んだ線分のベクトル

    wx2 = np.dot(w, x_vector) 
    xl2 = np.linalg.norm(x_vector)

    dx_array = wx2 / (wl * wl) * w # 正射影ベクトル
    dx = np.linalg.norm(dx_array)
    dy = np.sqrt((xl2 * xl2) - (dx * dx))
    
    if wx2/ (wl * wl) < 0 :
        dx = dx * -1

    if wx_cross < 0: # 外積の正負により角の向きを判定
        fi_deg = fi_deg * -1
        dy = dy * -1
    
    print(dx)
    print(dy)
    print(fi_deg)
    print("------")

offset(5, 2, 20, 10, 10, 10, 3.1415 / 2)

offset(5, 2, 20, 10, 25, 10, 0)




