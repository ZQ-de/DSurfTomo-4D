import os
import math
import numpy as np


# ============================================================
# 输入输出文件
# ============================================================
INPUT_FILE = "surfdata.dat"
OUTPUT_T1 = "surfdata_T1.dat"
OUTPUT_T2 = "surfdata_T2.dat"


# ============================================================
# 低速体范围
# ============================================================
LAT_MIN = 25.08
LAT_MAX = 25.12
LON_MIN = 121.43
LON_MAX = 121.48

# 低速体内部速度系数
# 0.90 = 低速体内部速度降低 10%
VELOCITY_FACTOR = 0.90


# ============================================================
# 局部经纬度转平面坐标
# 小范围内用等距近似即可
# x: km, y: km
# ============================================================
LAT_REF = (LAT_MIN + LAT_MAX) / 2.0
LON_REF = (LON_MIN + LON_MAX) / 2.0


def lonlat_to_xy_km(lon, lat):
    """
    Convert lon/lat to local Cartesian coordinates in km.
    """
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * math.cos(math.radians(LAT_REF))

    x = (lon - LON_REF) * km_per_deg_lon
    y = (lat - LAT_REF) * km_per_deg_lat

    return x, y


def segment_box_intersection_length(x1, y1, x2, y2, xmin, xmax, ymin, ymax):
    """
    Liang-Barsky clipping.
    Return the length of the segment inside the rectangle.
    """

    dx = x2 - x1
    dy = y2 - y1

    p = [-dx, dx, -dy, dy]
    q = [x1 - xmin, xmax - x1, y1 - ymin, ymax - y1]

    u1 = 0.0
    u2 = 1.0

    for pi, qi in zip(p, q):
        if abs(pi) < 1e-12:
            if qi < 0:
                return 0.0
        else:
            r = qi / pi
            if pi < 0:
                if r > u2:
                    return 0.0
                if r > u1:
                    u1 = r
            else:
                if r < u1:
                    return 0.0
                if r < u2:
                    u2 = r

    if u2 <= u1:
        return 0.0

    total_len = math.sqrt(dx * dx + dy * dy)
    inside_len = (u2 - u1) * total_len

    return inside_len


def generate_4d_lowbody():
    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(f"找不到输入文件: {INPUT_FILE}")

    # 低速体矩形边界转平面坐标
    x_min, y_min = lonlat_to_xy_km(LON_MIN, LAT_MIN)
    x_max, y_max = lonlat_to_xy_km(LON_MAX, LAT_MAX)

    xmin, xmax = min(x_min, x_max), max(x_min, x_max)
    ymin, ymax = min(y_min, y_max), max(y_min, y_max)

    total_rays = 0
    affected_rays = 0

    inside_lengths = []
    path_lengths = []
    dv_list = []
    dt_list = []
    ratio_list = []

    src_lat = None
    src_lon = None

    with open(INPUT_FILE, "r") as fin, \
         open(OUTPUT_T1, "w") as f1, \
         open(OUTPUT_T2, "w") as f2:

        for raw_line in fin:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith("#"):
                f1.write(line + "\n")
                f2.write(line + "\n")

                parts = line.split()
                src_lat = float(parts[1])
                src_lon = float(parts[2])

                continue

            parts = line.split()

            rec_lat = float(parts[0])
            rec_lon = float(parts[1])
            v1 = float(parts[2])

            x1, y1 = lonlat_to_xy_km(src_lon, src_lat)
            x2, y2 = lonlat_to_xy_km(rec_lon, rec_lat)

            dx = x2 - x1
            dy = y2 - y1
            path_len = math.sqrt(dx * dx + dy * dy)

            if path_len <= 0 or v1 <= 0:
                v2 = v1
                inside_len = 0.0
                dt = 0.0
            else:
                inside_len = segment_box_intersection_length(
                    x1, y1, x2, y2,
                    xmin, xmax, ymin, ymax
                )

                if inside_len > 1e-8:
                    # 原始走时
                    t1 = path_len / v1

                    # 低速体内部速度
                    v_inside = v1 * VELOCITY_FACTOR

                    # 只让穿过低速体的那一段变慢
                    t2 = (path_len - inside_len) / v1 + inside_len / v_inside

                    # 再转回 DSurfTomo 需要的等效相速度
                    v2 = path_len / t2

                    dt = t2 - t1
                    affected_rays += 1

                    inside_lengths.append(inside_len)
                    path_lengths.append(path_len)
                    dv_list.append(v2 - v1)
                    dt_list.append(dt)
                    ratio_list.append(inside_len / path_len)
                else:
                    v2 = v1
                    dt = 0.0

            total_rays += 1

            f1.write(f" {rec_lat:10.6f} {rec_lon:10.6f} {v1:10.4f}\n")
            f2.write(f" {rec_lat:10.6f} {rec_lon:10.6f} {v2:10.4f}\n")

    print("✅ 低速体 T1/T2 数据生成完成")
    print(f"输入文件: {INPUT_FILE}")
    print(f"输出 T1: {OUTPUT_T1}")
    print(f"输出 T2: {OUTPUT_T2}")

    print("\n低速体范围:")
    print(f"lat: {LAT_MIN} – {LAT_MAX}")
    print(f"lon: {LON_MIN} – {LON_MAX}")
    print(f"velocity_factor = {VELOCITY_FACTOR}")
    print(f"低速体内部速度降低比例 = {(1.0 - VELOCITY_FACTOR) * 100:.1f}%")

    print("\n射线统计:")
    print(f"总射线数: {total_rays}")
    print(f"穿过低速体的射线数: {affected_rays}")
    print(f"未受影响射线数: {total_rays - affected_rays}")

    if affected_rays > 0:
        inside_lengths = np.array(inside_lengths)
        path_lengths = np.array(path_lengths)
        dv_list = np.array(dv_list)
        dt_list = np.array(dt_list)
        ratio_list = np.array(ratio_list)

        print("\n穿过低速体的路径长度统计，单位 km:")
        print(f"Lin min  = {inside_lengths.min():.4f}")
        print(f"Lin max  = {inside_lengths.max():.4f}")
        print(f"Lin mean = {inside_lengths.mean():.4f}")
        print(f"Lin std  = {inside_lengths.std():.4f}")

        print("\n低速体路径占比 Lin/L 统计:")
        print(f"min  = {ratio_list.min():.4f}")
        print(f"max  = {ratio_list.max():.4f}")
        print(f"mean = {ratio_list.mean():.4f}")
        print(f"std  = {ratio_list.std():.4f}")

        print("\n速度差 V2 - V1 统计，单位 km/s:")
        print(f"min  = {dv_list.min():.6f}")
        print(f"max  = {dv_list.max():.6f}")
        print(f"mean = {dv_list.mean():.6f}")
        print(f"std  = {dv_list.std():.6f}")

        print("\n走时差 T2 - T1 统计，单位 s:")
        print(f"min  = {dt_list.min():.6f}")
        print(f"max  = {dt_list.max():.6f}")
        print(f"mean = {dt_list.mean():.6f}")
        print(f"std  = {dt_list.std():.6f}")
        print(f"rms  = {np.sqrt(np.mean(dt_list ** 2)):.6f}")

    else:
        print("\n⚠️ 没有射线穿过低速体。请调整低速体范围。")


if __name__ == "__main__":
    generate_4d_lowbody()
