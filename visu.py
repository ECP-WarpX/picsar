import numpy as np
import matplotlib.pyplot as plt


w0=187106056556544.47
dt=2.324281450142292e-16
em=np.loadtxt("em.txt")
ek=np.loadtxt("kin_e.txt")
plt.ion()
plt.title("em")
plt.plot(em)
plt.figure()
plt.plot(ek)
plt.title("ek")
et=ek+em
plt.figure()
plt.title("et")
plt.plot(et)
