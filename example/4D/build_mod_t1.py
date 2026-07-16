import numpy as np


T1_MODEL_FILE = "Result_T1/DSurfTomo.inMeasure.dat"
MOD_TEMPLATE_FILE = "MODini"
OUT_MOD_FILE = "MOD_T1"

NX = 18
NY = 18
NZ = 9

LAT0 = 25.200
LON0 = 121.350
DLAT = 0.015
DLON = 0.017


def read_mod_template(path):
    with open(path, "r") as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    depths = np.array([float(x) for x in lines[0].split()])
    data_lines = lines[1:]

    expected = NZ * NY
    if len(data_lines) != expected:
        raise ValueError(f"MOD 数据行数不对：应为 {expected}，实际 {len(data_lines)}")

    # 对应 Fortran: vsf(i,j,k)
    # 文件顺序: k -> j -> i
    # Python 保存为 model[k, j, i]
    model = np.zeros((NZ, NY, NX), dtype=float)

    idx = 0
    for k in range(NZ):
        for j in range(NY):
            vals = [float(x) for x in data_lines[idx].split()]
            if len(vals) != NX:
                raise ValueError(f"第 {idx + 2} 行列数不对：应为 {NX}，实际 {len(vals)}")
            model[k, j, :] = vals
            idx += 1

    return depths, model


def write_mod(path, depths, model):
    with open(path, "w") as f:
        f.write("  " + "  ".join([f"{z:.1f}" for z in depths]) + "\n")

        for k in range(NZ):
            for j in range(NY):
                line = "".join([f"{model[k, j, i]:7.3f}" for i in range(NX)])
                f.write(line + "\n")


def build_mod_t1():
    print("读取 MODini 模板...")
    depths, model = read_mod_template(MOD_TEMPLATE_FILE)

    print("读取 T1 最终模型...")
    t1 = np.loadtxt(T1_MODEL_FILE)
    if t1.ndim == 1:
        t1 = t1.reshape(1, -1)

    replaced = 0
    skipped = 0

    for row in t1:
        lon, lat, dep, vs = row[:4]

        # DSurfTomo 输出:
        # lon = LON0 + (j_out - 1) * DLON
        # lat = LAT0 - (i_out - 1) * DLAT
        #
        # 但输出速度是:
        # vsf(i_out + 1, j_out + 1, k)
        #
        # 所以写回 MOD 时必须整体向内偏移一格。
        i_out = int(round((LAT0 - lat) / DLAT))
        j_out = int(round((lon - LON0) / DLON))

        i_mod = i_out + 1
        j_mod = j_out + 1

        k_mod = int(np.argmin(np.abs(depths - dep)))

        if not (0 <= i_mod < NX and 0 <= j_mod < NY and 0 <= k_mod < NZ):
            skipped += 1
            continue

        if abs(depths[k_mod] - dep) > 1e-4:
            skipped += 1
            continue

        model[k_mod, j_mod, i_mod] = vs
        replaced += 1

    write_mod(OUT_MOD_FILE, depths, model)

    print("✅ MOD_T1 生成完成")
    print("T1 输入节点数:", len(t1))
    print("成功替换节点数:", replaced)
    print("跳过节点数:", skipped)
    print("输出文件:", OUT_MOD_FILE)
    print("深度层:", depths)

    with open(OUT_MOD_FILE, "r") as f:
        nline = sum(1 for _ in f)

    print("MOD_T1 总行数:", nline)
    print("理论总行数:", 1 + NZ * NY)

    if nline != 1 + NZ * NY:
        raise ValueError("MOD_T1 行数错误")

    print("✅ MOD_T1 行数检查通过")


if __name__ == "__main__":
    build_mod_t1()
