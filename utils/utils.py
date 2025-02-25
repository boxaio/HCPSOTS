import numpy as np
import torch
import random
import matplotlib.pyplot as plt






def compute_fft(x):
    x_0 = torch.fft.fftshift(torch.fft.fft2(x[:, 0, :, :]))
    x_1 = torch.fft.fftshift(torch.fft.fft2(x[:, 1, :, :]))
    x_2 = torch.fft.fftshift(torch.fft.fft2(x[:, 2, :, :]))
    x = torch.stack([x_0, x_1, x_2], dim=1)
    return x





batch_size = 1
device = torch.device("cuda:0")
white_noise = torch.randn(batch_size, 3, 256, 256).float().to(device)
white_noise_plot = (white_noise - white_noise.min()) / (white_noise.max() - white_noise.min())

white_noise_fft = torch.abs(compute_fft(white_noise))

blue_pts = np.load('/media/box/Elements/Exp/HCPSOTS/out/blue_points_w500_h500.npz')
blue_pts = blue_pts['blue_points'].flatten()
blue_pts = 2*(blue_pts - blue_pts.min()) / (blue_pts.max() - blue_pts.min())-1

# print(blue_pts.min(), blue_pts.max())

nrange = list(np.arange(blue_pts.shape[0]))

gaussian_blue_noise = blue_pts[random.sample(nrange, batch_size*3*64*64)]
print(gaussian_blue_noise.mean(), gaussian_blue_noise.std())
gaussian_blue_noise = torch.tensor(gaussian_blue_noise.reshape(batch_size, 3, 64, 64)).float().to(device)

gaussian_blue_noise_plot = (gaussian_blue_noise - gaussian_blue_noise.min()) \
                              / (gaussian_blue_noise.max() - gaussian_blue_noise.min())
gaussian_blue_noise_fft = torch.abs(compute_fft(gaussian_blue_noise))




output_dir = "/media/box/Elements/Exp/HCPSOTS/out"

plt.figure(figsize=(20,18))
plt.subplot(221)
plt.title('white noise')
plt.imshow(white_noise_plot[0].permute(1, 2, 0).cpu().numpy())
plt.subplot(222)
plt.title('white noise spectrum')
plt.imshow(white_noise_fft[0, 0].real.detach().cpu().numpy(), cmap='gray')
plt.subplot(223)
plt.title('blue noise')
plt.imshow(gaussian_blue_noise_plot[0].permute(1, 2, 0).cpu().numpy())
plt.subplot(224)
plt.title('blue noise spectrum')
plt.imshow(gaussian_blue_noise_fft[0, 0].real.detach().cpu().numpy(), cmap='gray')
plt.tight_layout()
plt.savefig(output_dir+'/noise_spectrum.png', dpi=600)
plt.clf()

