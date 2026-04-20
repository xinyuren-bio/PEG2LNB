#!/usr/bin/env python3
"""
消除 GRO 文件中严重重叠的原子对，使最小化能正常起步。
仅处理距离 < clash_nm 的原子对，将二者沿连线方向推开至 target_nm。
考虑周期边界条件 (PBC)。
不依赖 numpy，仅用标准库。
"""
import sys
import math

def read_gro(path):
    with open(path) as f:
        title = f.readline().strip()
        n = int(f.readline().strip())
        lines = []
        for _ in range(n):
            lines.append(f.readline())
        box_line = f.readline()
    coords = []
    for line in lines:
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        coords.append([x, y, z])
    box = [float(x) for x in box_line.split()[:3]]
    return title, n, lines, coords, box

def write_gro(path, title, n, lines, coords, box):
    with open(path, 'w') as f:
        f.write(title + '\n')
        f.write('%d\n' % n)
        for i in range(n):
            x, y, z = coords[i]
            new_line = lines[i][:20] + '%8.3f%8.3f%8.3f' % (x, y, z) + lines[i][44:]
            f.write(new_line)
        f.write('%10.5f%10.5f%10.5f' % tuple(box[:3]))
        if len(box) > 3:
            f.write(''.join(['%10.5f' % b for b in box[3:]]))
        f.write('\n')

def pbc_dx(dx, b):
    return dx - round(dx / b) * b

def pbc_vector(ri, rj, box):
    dx = pbc_dx(rj[0] - ri[0], box[0])
    dy = pbc_dx(rj[1] - ri[1], box[1])
    dz = pbc_dx(rj[2] - ri[2], box[2])
    return [dx, dy, dz]

def norm(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def fix_overlaps(coords, box, clash_nm=0.20, target_nm=0.28):
    n = len(coords)
    total_moves = 0
    # 邻格搜索 + 按格子字典序去重：每对原子在整个 round 内只检查一次（避免旧实现中同一对被重复处理数万次）
    # 膜/界面附近局部密度高，格子过粗会导致单个邻域内原子过多；1 nm 量级较稳妥
    cell_size = max(1.0, min(box) / 25.0)
    nx = max(1, int(box[0] / cell_size))
    ny = max(1, int(box[1] / cell_size))
    nz = max(1, int(box[2] / cell_size))
    cells = {}
    cell_of = []
    for i in range(n):
        c = coords[i]
        ix = int(c[0] / box[0] * nx) % nx
        iy = int(c[1] / box[1] * ny) % ny
        iz = int(c[2] / box[2] * nz) % nz
        key = (ix, iy, iz)
        if key not in cells:
            cells[key] = []
        cells[key].append(i)
        cell_of.append((ix, iy, iz))
    for i in range(n):
        ix, iy, iz = cell_of[i]
        ci = (ix, iy, iz)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    cx = (ix + di) % nx
                    cy = (iy + dj) % ny
                    cz = (iz + dk) % nz
                    ck = (cx, cy, cz)
                    if ck not in cells:
                        continue
                    if ck > ci:
                        js = cells[ck]
                    elif ck == ci:
                        js = [j for j in cells[ck] if j > i]
                    else:
                        continue
                    for j in js:
                        v = pbc_vector(coords[i], coords[j], box)
                        d = norm(v)
                        if d < 1e-6:
                            d = 1e-6
                        if d < clash_nm:
                            move = (target_nm - d) / 2.0
                            u = [v[0]/d, v[1]/d, v[2]/d]
                            for k in range(3):
                                coords[i][k] -= move * u[k]
                                coords[j][k] += move * u[k]
                            for k in range(3):
                                coords[i][k] = coords[i][k] % box[k]
                                coords[j][k] = coords[j][k] % box[k]
                            total_moves += 1
    return total_moves

def main():
    clash_nm = 0.35
    target_nm = 0.48
    max_rounds = 15
    if len(sys.argv) < 2:
        print('Usage: python fix_overlap_gro.py <input.gro> [output.gro]')
        print('  Default output: input_fixed.gro')
        sys.exit(1)
    inp = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else inp.replace('.gro', '_fixed.gro')
    if out == inp:
        out = inp.replace('.gro', '_fixed.gro')

    title, n, lines, coords, box = read_gro(inp)
    print('Atoms: %d, box: %s' % (n, box[:3]))
    for r in range(max_rounds):
        moves = fix_overlaps(coords, box, clash_nm=clash_nm, target_nm=target_nm)
        if moves == 0:
            print('Round %d: no pairs < %.2f nm, done.' % (r + 1, clash_nm))
            break
        print('Round %d: resolved %d overlapping pairs' % (r + 1, moves))
    write_gro(out, title, n, lines, coords, box)
    print('Written: %s' % out)

if __name__ == '__main__':
    main()
