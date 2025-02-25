import numpy as np
import matplotlib.pyplot as plt
from noise import PoissonDisc, UniformNoise

# Generate periodogram images for uniformly-distributed noise and
# Poisson disc-sampled ("blue") noise in two dimensions.
# For mathematical details,  please see the blog articles at
# https://scipython.com/blog/poisson-disc-sampling-in-python/
# https://scipython.com/blog/power-spectra-for-blue-and-uniform-noise/
# Christian Hill, March 2017.



# domain size, minimum distance between samples for Poisson disc method...
width = height = 100
r = 2
poisson_disc = PoissonDisc(width, height, r)
# Expected number of samples from Poisson disc method...
n = int(width * height / np.pi / poisson_disc.a**2)

# ... use the same for uniform noise.
uniform_noise = UniformNoise(width, height, n)

# Number of sampling runs to do (to remove noise from the noise in the power spectrum).
N = 10
# Sampling parameter, when putting the sample points onto the domain
M = 3

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16, 16))

for j, noise in enumerate((poisson_disc, uniform_noise)):
    print(noise.__class__.__name__)
    spectrum = np.zeros((height * M, width * M))
    for i in range(N):
        noise.reset()
        samples = np.array(noise.sample())
        diffs = np.diff(np.sort(samples[:,1]))
        print('{}/{}'.format(i+1, N), samples[:,0].min(), samples[:,0].max(), 
              np.abs(diffs).mean(), np.abs(diffs).var())
        domain = np.zeros((height * M, width * M))
        for pt in samples:
            coords = int(pt[1] * M), int(pt[0] * M)
            domain[coords] = 1.0

        # 2d Fourier Trasform, shift the frequencies 
        f = np.fft.fft2(domain)
        fshift = np.fft.fftshift(f)
        spectrum += np.log(np.abs(fshift))
    
    # Plot the a set of random points and the power spectrum.
    ax[0][j].imshow(domain, cmap=plt.cm.Greys, vmin=0, vmax=1)
    # ax[0][j].imshow(domain, extent=[0, width, 0, height])
    ax[1][j].imshow(spectrum, cmap=plt.cm.Greys_r)
    # Remove axis ticks and annotations
    for k in (0,1):
        ax[k][j].tick_params(which='both', bottom='off', left='off',
                top='off', right='off', labelbottom='off', labelleft='off')

plt.tight_layout()
plt.savefig('/media/box/Elements/Exp/HCPSOTS/out/periodograms.png', dpi=600)
# plt.show()
