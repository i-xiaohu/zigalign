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
        s = f.readline().split()
        len1 = int(s[2].rstrip(','))
        len2 = int(s[5])
        f.readline()
        bp1 = [int(x) for x in f.readline().split()]
        f.readline()
        bp2 = [int(x) for x in f.readline().split()]
        for line in f:
            x, y, c = line.strip().split()
            x, y, c = int(x), int(y), int(c)
            points.append((x, y, c))
    width = height = 8
    plt.figure(figsize=(width, height), dpi=350)
    ax = plt.gca()
    ax.set_title('Sequence Alignment With Duplications\n')
    ax.set_ylabel('Seq1 length: %d bp' % len1)
    ax.set_xlim(0, max([y for x, y, c in points]))
    ax.set_xlabel('Seq2 length: %d bp' % len2)
    ax.set_ylim(0, max([x for x, y, c in points]))
    ax.set_aspect('equal', adjustable='box')
    ax.invert_yaxis()
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    # ax.grid(True, color='gray', linestyle='-', alpha=0.3)

    bp_min1, bp_max1 = min(bp1), max(bp1)
    bp_min2, bp_max2 = min(bp2), max(bp2)
    for bp in bp1:
        plt.plot([bp_min2, bp_max2], [bp, bp], color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    for bp in bp2:
        plt.plot([bp, bp], [bp_min1, bp_max1], color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)

    lx, ly, lc = points[0]
    for (x, y, c) in points[1:]:
        # FIXME: just a guess; ZigAlign should provide accurate result by DP
        if lc == 3:
            if ly == y:
                offset = lx - x
                plt.plot([y-offset, y], [x, lx], color='blue', linestyle='--')
                plt.plot([y-offset, y], [x, x], color='blue', linestyle='--')
        if lc == 3:
            if ly == y:
                offset = lx - x
                plt.text(y, (lx + x) / 2, 'copy %d' % offset, fontsize=7)
            else:
                offset = ly - y
                plt.text(y, x, 'delete %d' % offset, fontsize=7)
                plt.plot([ly, y], [lx, x], color=color_list[lc])
        else:
            plt.plot([ly, y], [lx, x], color=color_list[lc])
        lx, ly, lc = x, y, c

    plt.plot([0, 0], [0, 0], color='blue', label='indel')
    plt.plot([0, 0], [0, 0], color='green', label='mat/mis')
    plt.plot([0, 0], [0, 0], color='grey', linestyle='--', label='breakpoint')
    plt.plot([0, 0], [0, 0], color='red', label='dup del')
    plt.plot([0, 0], [0, 0], color='blue', linestyle='--', label='copy')
    plt.legend(loc='upper right')

    fig_fn = bt_fn.rstrip(".txt") + ".png"
    plt.savefig(fig_fn)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: python vis_pair.py <prefix.txt>\n')
        sys.stderr.write('Visualized figure is stored as prefix.png\n')
    else:
        visualize(sys.argv[1])
