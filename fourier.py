"""
用于解析path并生成GIF图像
"""

import math
import cmath
import re
import time
import multiprocessing
from turtle import color
import matplotlib.pyplot as plt
from fourier.clean_data import get_data_from_func, clean_data
from fourier.draw import draw
from matplotlib.animation import FuncAnimation
import numpy as np

def bezier(t,a,b,c,d):#贝塞尔函数（估测长度用）
    return (-a+3*b-3*c+d)*t*t*t+3*(a-2*b+c)*t*t+3*(-a+b)*t+a

def linear(x,a,b,c,d):#映射函数
    return (x-a)/(b-a)*(d-c)+c
'''
def equation(bt,t,m,bV):#主要计算函数——（废弃）
    return bezier(bt,bV[0],bV[1],bV[2],bV[3])*cmath.exp(-m*1j*t)
'''
oneOver2PI = 1/(math.pi*2)#常量

def prSolve(m,cs,ce,n):#主要计算函数1-解析解
    if m == 0 :
        return (ce-cs)*oneOver2PI/(n+1)
    if n==0 :
        return 1j*oneOver2PI/m*(cmath.exp(-m*1j*ce)-cmath.exp(-m*1j*cs))
    elif n>0 :
        return 1j*oneOver2PI/m*cmath.exp(-m*1j*ce)-n*1j/((ce-cs)*m)*prSolve(m,cs,ce,n-1)
    else : 
        return 0

def numSolve(m,cs,ce,pts):#主要计算函数2-贝塞尔曲线方程代入
    return (-pts[0]+3*pts[1]-3*pts[2]+pts[3])*prSolve(m,cs,ce,3)+3*(pts[0]-2*pts[1]+pts[2])*prSolve(m,cs,ce,2)+3*(-pts[0]+pts[1])*prSolve(m,cs,ce,1)+pts[0]*prSolve(m,cs,ce,0)


def cpToList(cp):
    return [cp.real,cp.imag]
start = 0
end = 1000#圆圈数量
prjNum = end-start+1
prjCur = None
points = []#贝塞尔采集点
out = []#输出用坐标容器
center = [500,500]#中心点位置
curWeight = []#各段曲线的时间权重容器
trDatas = []#临时容器
sharedPoints = None

def initParam(pts,cw,pjc):
    global points,curWeight,prjCur
    points = pts
    curWeight = cw
    prjCur = pjc

def mainCalculation(s):
    global points,curWeight,prjCur
    m = 0
    if s>0: 
        m = ((s+1)//2)*(-1 if (s%2 == 0) else 1)
    #print("Now working on orbit {0},m = {1}".format(s,m))
    sum = 0j
    for i in range(len(points)):
        cs = linear(curWeight[i],0,1,0,math.pi*2)
        ce = linear(curWeight[i+1],0,1,0,math.pi*2)
        sum += numSolve(m,cs,ce,points[i])
    with prjCur.get_lock():
        prjCur.value += 1
    return cpToList(sum)


if __name__ == '__main__':
    with open('rawvertexes.txt', 'r') as f:#读取并解析SVG信息
        rawdata = f.readlines()
        curve = re.sub(r'\s','',"".join(rawdata))
        cells = re.findall(r'\w[\d\,\-\.]+',curve)
        for cell in cells:
            trcdata = []
            formatString = re.sub(r'-',',-',cell)
            trcdata.append(re.match(r'[A-Za-z]',formatString).group(0))
            rawvers = re.sub(r'[A-Za-z]\,?','',formatString).split(',')
            for st in range(len(rawvers)):
                rawvers[st] = float(rawvers[st])
            vergroup = []
            vercurgrp = []
            for st in range(len(rawvers)):
                vercurgrp.append(rawvers[st])
                if len(vercurgrp)>=2:
                    vergroup.append(vercurgrp[0]+vercurgrp[1]*1j)
                    vercurgrp.clear()
            if len(vercurgrp)>0:
                if re.match(r'v',trcdata[0],re.I):
                    vergroup.append(0+vercurgrp[0]*1j)
                elif re.match(r'h',trcdata[0],re.I):
                    vergroup.append(vercurgrp[0]+0j)
            trcdata.append(vergroup)
            trDatas.append(trcdata)
        print('Vertexes data have been read.')
        print(type(trcdata))
        print(trDatas)

    for i in range(1,len(trDatas)):#解析SVG信息
        if re.match(r'[a-z]',trDatas[i][0]):
            for j in range(len(trDatas[i][1])):
                trDatas[i][1][j]+=trDatas[i-1][1][-1]
    for i in range(len(trDatas)):#解析SVG信息
        flag = trDatas[i][0]
        if re.match(r'm',flag,re.I):continue
        trDatas[i][1].insert(0,trDatas[i-1][1][-1])
        if re.match(r's',flag,re.I):
            trDatas[i][1].insert(1,2*trDatas[i-1][1][-1]-trDatas[i-1][1][-2])
        if re.match(r'[lvh]',flag,re.I):
            trDatas[i][1].insert(1,trDatas[i][1][0]/3+trDatas[i][1][-1]*2/3)
            trDatas[i][1].insert(1,trDatas[i][1][0]*2/3+trDatas[i][1][-1]/3)

    for i in range(len(trDatas)):#解析SVG信息
        flag = trDatas[i][0]
        if re.match(r'm',flag,re.I):continue
        points.append(trDatas[i][1])
    for i in range(len(points)):#将中心点归零
        for j in range(len(points[i])):
            points[i][j] -= center[0]+1j*center[1]
    print(points)
    

    print("Weight process start.")#计算时间权重
    wsum = 0
    for curve in points:#Calculate weight
        wst = 10 #steps
        
        sum = 0
        for i in range(1,wst):
            sum+=abs(bezier(linear(i,0,wst,0,1),curve[0],curve[1],curve[2],curve[3])-
            bezier(linear(i-1,0,wst,0,1),curve[0],curve[1],curve[2],curve[3]))
        curWeight.append(sum)
        wsum+=sum
    for i in range(len(curWeight)):
        curWeight[i] /= wsum
    for i in range(1,len(curWeight)):
        curWeight[i] += curWeight[i-1]
    curWeight.insert(0,0)
    curWeight[-1] = 1
    print("Weight process finished.")
    print(points)
    print(len(points))
    print(points[0])
    print(type(points[1]))
    print(len(points[0]))
    print(type(points[1][0]))
    print(points[1][1])
    
    # 储存这些复数
    datas = []
    for i in range(len(points)):
        for j in range(len(points[0])):                                                                                                                                                                                                                                       
            datas.append(points[i][j])
    datas = np.array(datas)
    print(datas)
    print(type(datas)) 
    
    X = clean_data(datas)
    fig, ax = plt.subplots(1, 1,figsize=(20,10),dpi = 200)
    ax.set_facecolor('black')
    fig.set_facecolor('black') 
    update_all = draw(X, fig, ax)
    ani = FuncAnimation(fig, update_all, blit=True, interval=3, frames=len(X))
    # plt.show()
    ani.save('demo2.gif', writer='pillow')
    
    # ## 画最终的图
    # x = []
    # y = []
    # fig2, ax2 = plt.subplots(1, 1,figsize=(20,10))
    # for i in range(len(points)):
    #     for j in range(len(points[0])):
    #         x.append(points[i][j].real)
    #         y.append(points[i][j].imag)
    # ax2.set_xlim(-2600, 2600)
    # ax2.set_ylim(-1300, 1300)
    # ax2.set_xticks([])
    # ax2.set_yticks([]) 
    # ax2.set_facecolor('black')
    # fig2.set_facecolor('black')
    # ax2.plot(x,y,color="white",linestyle = '-')     
    # plt.savefig("last.jpg-")
    

