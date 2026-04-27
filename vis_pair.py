import sys
import matplotlib.pyplot as plt


# Input: backtrace txt file
# Output: figure png file
def visualize(bt_fn):
    points = []
    # colors annotation:
    # 0, blue: gap
    # 1, blue: gap
    # 2, green: match/mismatch
    # 3, red: duplication deletions
    color_list = ['blue', 'blue', 'green', 'red']
    with open(bt_fn, 'r') as f:
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

    lx, ly, lc = points[0]
    for (x, y, c) in points[1:]:
        plt.plot([ly, y], [lx, x], color=color_list[lc])
        lx, ly, lc = x, y, c

    plt.plot([0, 0], [0, 0], color='blue', label='indel')
    plt.plot([0, 0], [0, 0], color='green', label='mat/mis')
    plt.plot([0, 0], [0, 0], color='red', label='dup del')
    plt.legend()

    fig_fn = bt_fn.rstrip(".txt") + ".png"
    plt.savefig(fig_fn)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: python vis_pair.py <prefix.txt>\n')
        sys.stderr.write('Visualized figure is stored as prefix.png\n')
    else:
        visualize(sys.argv[1])
