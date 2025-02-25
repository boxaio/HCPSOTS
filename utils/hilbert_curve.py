import numpy as np
import matplotlib.pyplot as plt



def hilbert_curve(x0, y0, xi, yi, xj, yj, n):
    points = []
    if n <= 0:
        X = x0 + (xi + yi) / 2
        Y = y0 + (xj + yj) / 2
        points.append((X, Y))
    else:
        n -= 1
        points.extend(hilbert_curve(          x0,           y0,  yi/2,  xi/2,  yj/2,  xj/2, n))
        points.extend(hilbert_curve(     x0+xi/2,      y0+xj/2,  xi/2,  yi/2,  xj/2,  yj/2, n))
        points.extend(hilbert_curve(x0+xi/2+yi/2, y0+xj/2+yj/2,  xi/2,  yi/2,  xj/2,  yj/2, n))
        points.extend(hilbert_curve(  x0+xi/2+yi,   y0+xj/2+yj, -yi/2, -xi/2, -yj/2, -xj/2, n))
    return points


def curve_func(points: np.ndarray, num_pts_per_segment: int=5):
    """
    Define a continuous curve using the given Hilbert points.
    points: np.ndarray, shape=(N, 2), the Hilbert points
    num_pts_per_segment: int, number of points per segment
    """
    N = points.shape[0]
    # lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)

    # segs = np.linspace(0.0, 1.0, num_pts_per_segment)
    segs = np.arange(0, 1.0, 1.0/num_pts_per_segment)

    c_pts = []
    for i in range(1, N):
        v = points[i] - points[i-1]
        for k in segs:
            c_pts.append(points[i-1] + v * k)

    return np.array(c_pts)


points_1 = np.array(hilbert_curve(0, 0, 1, 0, 0, 1, n=1))
points_2 = np.array(hilbert_curve(0, 0, 1, 0, 0, 1, n=2))
points_3 = np.array(hilbert_curve(0, 0, 1, 0, 0, 1, n=3))
# print(points_1.shape, points_2.shape, points_3.shape)



if __name__ == '__main__':

    c_pts = curve_func(points_3)
    print(c_pts)


    plt.figure(figsize=(10, 8), )
    plt.plot(c_pts[:, 0], c_pts[:, 1], '#236B8E')
    plt.plot(c_pts[:, 0], c_pts[:, 1], '.', color='firebrick')
    plt.axis('equal')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()
