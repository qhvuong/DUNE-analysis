import numpy as np 
from matplotlib import pyplot as plt 

#N_values = [120, 150, 200, 250, 350, 600]
N_values = [50,70,100,150,220,320,500,800]
dm2_values = [0.0023, 0.01, 1.0000]

fig, axs = plt.subplots(2, 4, figsize=(16, 8))

for i, ax in enumerate(axs.flat):
    N = N_values[i]
    print(N)
    x = np.linspace(0, 1200, N)

    #plt.figure(figsize=(10, 6))

    for dm2 in dm2_values:
        delta = 1.27 * x * dm2
        prob = np.sin(delta)**2
        ax.plot(x,prob,alpha=0.8,label='${\\Delta}m^2$ = 'f'{dm2:.4f}''$eV^2$')

    ax.set_title(f'Nbins = {N}')
    ax.legend(loc='upper right')
    ax.set_xlabel('$L/E$ [km/MeV]')
    ax.set_ylabel('probability')
    #ax.savefig(f"prob_{N}.png")

plt.tight_layout()

fig.savefig('prob.png')
#plt.show()
plt.close()


