import sys

import matplotlib.pyplot as plt
import re

def parse_cigar(cigar):
    """解析CIGAR字符串，返回操作列表 [(op, length), ...]"""
    pattern = r'(\d+)([MIDX])'
    return [(op, int(length)) for length, op in re.findall(pattern, cigar)]

def cigar_to_coords(cigar):
    """
    根据CIGAR生成两条序列的坐标点列表
    返回: (seq1_coords, seq2_coords)
    """
    ops = parse_cigar(cigar)
    s1, s2 = 0, 0
    points = [(0, 0)]

    for op, length in ops:
        if op == 'M' or op == 'X':
            s1 += length
            s2 += length
            points.append((s1, s2))
        elif op == 'D':
            s1 += length
            points.append((s1, s2))
        elif op == 'I':
            s2 += length
            points.append((s1, s2))

    return points

def plot_cigar(cigar, output_file):
    """
    根据CIGAR绘制共线性图（保留坐标轴，去掉图例/标题/信息框）
    """
    points = cigar_to_coords(cigar)
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]

    fig, ax = plt.subplots(figsize=(10, 8))

    # 绘制路径线
    ax.plot(xs, ys, 'b-', linewidth=1.0, alpha=0.8)

    # 标记起点和终点（小点）
    ax.plot(xs[0], ys[0], 'go', markersize=4)
    ax.plot(xs[-1], ys[-1], 'ro', markersize=4)

    # 对角线参考线（浅灰色虚线）
    # max_len = max(xs[-1], ys[-1])
    # ax.plot([0, max_len], [0, max_len], 'k--', alpha=0.15, linewidth=0.5)

    # ===== 保留坐标轴相关，去掉其它标签 =====
    ax.set_xlabel('Sequence 1')      # X轴标签
    ax.set_ylabel('Sequence 2')      # Y轴标签
    ax.set_title('')                 # 去掉标题
    ax.grid(True, alpha=0.2)         # 保留浅网格（方便读坐标）

    # 坐标轴刻度保留（自动）
    # 如果刻度太密集，可以手动设置间隔，例如：
    # from matplotlib.ticker import MultipleLocator
    # ax.xaxis.set_major_locator(MultipleLocator(10000))
    # ax.yaxis.set_major_locator(MultipleLocator(10000))

    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"图片已保存为 {output_file}")
    print(f"序列1长度: {xs[-1]}, 序列2长度: {ys[-1]}")

# ===== 测试 =====
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print('cigar2plot cigar.txt output.png')
        exit(0)

    with open(sys.argv[1], 'r') as f:
        cigar = f.readline().strip()
    plot_cigar(cigar, sys.argv[2])