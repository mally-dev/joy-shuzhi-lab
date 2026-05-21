import numpy as np
def householder_QR(A):
    Q=np.eye(A.shape[0])  # 初始化Q为单位矩阵
    R=A.copy()  # 初始化R为A的副本
    n=A.shape[0]
    for k in range(n):
        x=R[k:, k]  # 提取当前列的子向量
        e=np.zeros_like(x)
        sign =1 if x[0] >= 0 else -1  # 确定符号以避免数值不稳定
        e[0]=np.linalg.norm(x)  # 计算e的第一个分量
        v=x+sign*e  # 计算Householder向量
        w=v/np.linalg.norm(v)  # 规范化Householder向量
        H=np.eye(n)  # 初始化Householder矩阵为单位矩阵
        H[k:,k:]=np.eye(len(x))-2*np.outer(w, w)  # 更新Householder矩阵
        R=H @ R  # 更新R
        Q=Q @ H.T  # 更新Q
    return Q, R
n=int(input("Enter the size of the matrix: "))
A=np.zeros((n,n))
for i in range(n):
    row=input("Enter row {} of the matrix (space-separated): ".format(i+1))
    A[i] = [float(x) for x in row.split()]
# Q, R = householder_QR(A)
# print("Q:\n", Q)
# print("R:\n", R)
#givens略
# def Givens_QR():
#     Q=np.eye(A.shape[0])  # 初始化Q为单位矩阵
#     n=A.shape[0]
#     for k in 