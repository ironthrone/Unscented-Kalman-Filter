import matplotlib.pyplot as plt
import numpy as np

lidar_line = 5.99
radar_line = 7.82


def plot_nis(f, threshold, outfile):
    nis_file = open(f)
    nis = []
    for i in nis_file:
        nis.append(float(i))

    plt.plot((0, len(nis)), (threshold, threshold))
    plt.plot(nis)
    # plt.show()
    plt.savefig(outfile)
    plt.clf()

plot_nis("lidar_nis_file.txt",lidar_line,"nis_lidar.png")
plot_nis("radar_nis_file.txt",radar_line,"nis_radar.png")
