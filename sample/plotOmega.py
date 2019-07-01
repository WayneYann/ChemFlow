# Mon May 13 10:25:04 CST 2019
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#Files
figname = '1'
print("Enter file names:")
filename = []
while True:
    name_input = input('> ')
    if name_input == '':
        break
    else:
        filename.append(name_input)

#Basic settings
figs = (14,9)
fonts1 = 20
linew = 2
plt.rcParams['font.family'] = 'serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'

#Plot data
file = filename[0]
with open(file) as f:
    line = f.readline()
    line = line.split(',')
    name = []
    for i in range(len(line)):
        name.append(line[i])

fig1 = plt.figure(figsize=figs)
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)
for i,file in enumerate(filename):
    data = np.loadtxt(file,delimiter=',',comments='#',skiprows=1)
    data = np.transpose(data)
    x = data[0]
    for k in range(len(name)):
        if k in [2,7,9,10,18,19,20]:
            ax1.plot(x,data[k],label=name[k],ls='-',lw=linew)
    ax2.plot(x,data[1],label=name[1],ls='-',c='k',lw=linew+1)


#Fig
ax1.legend(loc=0, fontsize=fonts1)
ax2.legend(loc=0, fontsize=fonts1)
x_ticks = np.linspace(x[0], x[-1], 5)
ax1.set_xticks(x_ticks)
ax1.tick_params(labelsize=fonts1)
ax1.set_xlabel(r'$x$ (m)',fontsize=fonts1,color='k')
ax1.set_ylabel(r'$\dot{\omega}$ (kg/m3 s)',fontsize=fonts1,color='k')

ax2.set_xticks(x_ticks)
ax2.tick_params(labelsize=fonts1)
ax2.set_xlabel(r'$x$ (m)',fontsize=fonts1,color='k')
ax2.set_ylabel(r'$\dot{Q}$ (J/m3 s)',fontsize=fonts1,color='k')
plt.savefig(figname+'.png',dpi=500,bbox_inches='tight')
plt.show()
