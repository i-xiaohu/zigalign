import sys
import matplotlib.pyplot as plt


def visualize(test_id: int):
    points = []
    color_list = ['green', 'red', 'red', 'blue', 'orange']
    with open('self_%d.txt' % test_id, 'r') as f:
        for line in f:
            x, y, c = line.strip().split()
            x, y, c = int(x), int(y), int(c)
            points.append((x, y, c))

    width = height = 8
    plt.figure(figsize=(width, height), dpi=350)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.invert_yaxis()
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.grid(True, color='gray', linestyle='-', alpha=0.3)
    ax.set_axisbelow(True)
    plt.title('Case %d' % (test_id-1), fontsize=14, pad=30)

    lx, ly, lc = points[0]
    for (x, y, c) in points[1:]:
        plt.plot([ly, y], [lx, x], color=color_list[lc])
        lx, ly, lc = x, y, c

    plt.plot([0, 0], [0, 0], color='green', label='non-rep')
    plt.plot([0, 0], [0, 0], color='red', label='D transfer')
    plt.plot([0, 0], [0, 0], color='blue', label='copy')
    plt.plot([0, 0], [0, 0], color='orange', label='B transfer')
    plt.legend()

    plt.savefig('self_%d.png' % test_id)


if __name__ == '__main__':
    visualize(int(sys.argv[1]))
