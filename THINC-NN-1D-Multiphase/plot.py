import numpy as np
import matplotlib.pyplot as plt
import os

def plot(filename):
    data = np.loadtxt(filename, delimiter=',', skiprows=1)


    plt.figure(1)

    title = filename.rstrip('.csv')
    title = title.replace('-',' ')

    plt.suptitle(title, color='red')

    plt.subplot(2, 2, 1)
    plt.plot(data[:,0],data[:,1])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\rho$')
    plt.title('Density')
    plt.grid()

    plt.subplot(2, 2, 2)
    plt.plot(data[:,0],data[:,2],color='orange')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u$')
    plt.title('Velocity')
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.plot(data[:,0],data[:,3], color='green')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$p$')
    plt.title('Pressure')
    plt.grid()

    plt.subplot(2, 2, 4)
    plt.plot(data[:,0],data[:,4],color='purple')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\alpha_1$')
    plt.title('Volume Fraction')
    plt.grid()

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    
    outfile = filename.replace('csv','png')
    plt.savefig(outfile, dpi = 1024)

    plt.clf()


    return 0


for in_file in sorted(os.listdir('.')):
    if in_file.endswith('.csv'):
        print(in_file)
        plot(in_file)
