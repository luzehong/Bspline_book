import matplotlib.pyplot as plt #插值模块
import numpy as np #数值计算模块
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interpolate

file = 'cmake-build-debug/01.txt'
a = np.loadtxt(file)
x = a[:,0]
y = a[:,1]
z = a[:,2]

#for i,j in zip(x,y):
#    print(i,j)
fig=plt.figure()
ax=Axes3D(fig)
ax.scatter(x,y,z,c='r')


#plt.plot(x,y,'-',color='red')     #定义样式
#plt.xlabel('Wavelength(μm)')  #横坐标含义
#plt.ylabel('Data Value')  #纵坐标含义
#plt.title('Kaolinite')
plt.show()

