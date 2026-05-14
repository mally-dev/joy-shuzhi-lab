# import numpy as np
# def power_method(A, x0, tol=1e-6, max_iter=1000):
#     z_old = x0 / np.linalg.norm(x0)  # 初始规格化
    
#     for k in range(max_iter):  # 用for循环控制迭代次数
#         y = A @ z_old
#         # 取绝对值最大的分量，保留符号
#         idx = np.argmax(np.abs(y))
#         lambda1 = y[idx]
        
#         z_new = y / lambda1  # 规格化
        
#         # 判断收敛：相邻两次z的差是否足够小
#         error = np.linalg.norm(z_new - z_old, np.inf)
        
#         if error < tol:
#             print(f"收敛于第{k+1}次迭代")
#             return z_new, lambda1
        
#         z_old = z_new  # 更新
    
#     print(f"达到最大迭代次数{max_iter}")
#     return z_new, lambda1
     
# # n=int(input("Enter the size of the matrix: "))
# # A=np.zeros((n,n))
# # for i in range(n):
# #     row=input("Enter row {} of the matrix (space-separated): ".format(i+1))
# #     A[i] = [float(x) for x in row.split()]
# # x0=np.random.rand(n)
# # eigenvector, eigenvalue=power_method(A, x0)
# # print(eigenvector,eigenvalue)
# # 可以加上这段验证代码
# A = np.array([[2,4,6],[3,9,15],[4,16,36]])
# x = np.array([0.18586754, 0.44603237, 1.0])
# lam = 43.87999715

# Ax = A @ x
# lam_x = lam * x

# print("A·x =", Ax)
# print("λ·x =", lam_x)
# print("误差 =", np.linalg.norm(Ax - lam_x))
def fabonaci(n):
    a=[1,1]
    for i in range(2,n):
        a.append(a[i-1]+a[i-2])
    print(a[-1])
fabonaci(1000)