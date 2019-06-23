# Mon May 13 10:25:04 CST 2019
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#Files
# figname = input('Figure name:\n')
figname = '1'
# print('File names:')
filename = []
filename = ['output.csv']
# while True:
#     name_input = input('> ')
#     if name_input == '':
#         break
#     else:
#         filename.append(name_input)
# xLabel = input('x label:\n')
# yLabel = input('y label:\n')

#Basic settings
figs = (10,10)
fonts1 = 18
linew = 2
plt.rcParams['font.family'] = 'serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'
#plt.style.use('ggplot')
#plt.grid(True, axis = 'y', color='k',ls='--',linewidth=0.6)
cs = ['k','r','c','m','b']
lstyl = ['--','-','-','-','-','-','-']
mkstyl = ['.','^','s','d','h']


#Plot data
file = filename[0]
with open(file) as f:
    line = f.readline()
    line = line.split(',')
    name = []
    for i in range(len(line)):
        name.append(line[i])
xLabel = name[0]
uLabel = name[1]
VLabel = name[2]
TLabel = name[3]
YLabel = name[4]

fig1 = plt.figure(figsize=figs)
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)
for i,file in enumerate(filename):
    data = np.loadtxt(file,delimiter=',',comments='#',skiprows=1)
    data = np.transpose(data)
    x = data[0]
    u = data[1]
    V = data[2]
    T = data[3]
    Y = data[4]
    ax1.plot(x,u,label='u',color=cs[i],ls='-',lw=linew)
    ax2.plot(x,V,label='V',color=cs[i],ls='-',lw=linew)
    # ax1.plot(x,y,label=file,ls=lstyl[i],lw=linew,color=cs[i])


#Fig
# ax1.legend(loc=0)
# ax2.legend(loc=0)
x_ticks = np.linspace(x[0], x[-1], 5)
ax1.set_xticks(x_ticks)
ax1.tick_params(labelsize=fonts1)
ax1.set_xlabel(rf'{xLabel}',fontsize=fonts1,color='k')
ax1.set_ylabel(rf'{uLabel}',fontsize=fonts1,color='k')

ax2.set_xticks(x_ticks)
ax2.tick_params(labelsize=fonts1)
ax2.set_xlabel(rf'{xLabel}',fontsize=fonts1,color='k')
ax2.set_ylabel(rf'{VLabel}',fontsize=fonts1,color='k')

plt.savefig(figname+'.png',dpi=500,bbox_inches='tight')
plt.show()
