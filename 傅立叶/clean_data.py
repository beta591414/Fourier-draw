"""
处理数据，将其变为复数对
"""
import numpy as np


def get_data_from_pic():
    pass


def get_data_from_csv():
    pass


def get_data_from_func():
    n = 100
    t = np.linspace(0, 2 * np.pi, n)
    x = 16 * (np.sin(t)) ** 3
    y = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
    X = x + 1j * y*2
    return X


def sort_data(X):
    '''
    假设输入数据已经是有序的了
    '''
    return X


def clean_data(X):
    X = sort_data(X)
    x, y = X.real, -X.imag
    x = (x - x.min()) / (x.max() - x.min()) - 0.5
    y = (y - y.min()) / (y.max() - y.min()) - 0.5
    X = x + 1j * y*1
    return X
