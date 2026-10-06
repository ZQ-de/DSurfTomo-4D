#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
把单期T1反演结果 DSurfTomo.inMeasure.dat 转换为 DSurfTomo-4D 所需的 MOD。

默认在当前 pds1 目录运行：
    python3 make_4d_MOD_from_T1.py

默认输入：
    DSurfTomo.inMeasure.dat   # T1第10次迭代最终绝对速度模型，通常16×16×8
    MOD                       # 单期反演所用初始MOD，结构模板18×18×9

默认输出：
    MOD_4D                    # 供4D目录使用，复制后改名为MOD
    MOD_4D_conversion.txt     # 转换检查记录

转换原则：
    1. 从 DSurfTomo.inMeasure.dat 读取内部有效网格：
       longitude latitude depth Vs
    2. 将16×16×8反演结果写入18×18×9 MOD模板的内部节点；
    3. 保留MOD模板最外一圈边界节点；
    4. 保留反演输出中没有的2.5 km底层；
    5. 输出格式严格为：
       第一行：9个深度节点
       后续：按 depth(k) → longitude(j) → latitude(i) 排列，每行18个Vs。
"""

from pathlib import Path
import hashlib
import numpy as np


INPUT_MODEL = Path("DSurfTomo.inMeasure.dat")
TEMPLATE_MOD = Path("MOD")
OUTPUT_MOD = Path("MOD_4D")
REPORT_FILE = Path("MOD_4D_conversion.txt")


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_file(path: Path, description: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"找不到{description}：{path.resolve()}")
    return path


def read_mod_template(path: Path):
    lines = require_file(path, "MOD模板").read_text(
        encoding="utf-8",
        errors="ignore",
    ).splitlines()

    if len(lines) < 2:
        raise ValueError("MOD文件内容不足。")

    depths = np.asarray(
        [float(value) for value in lines[0].split()],
        dtype=float,
    )

    values = np.asarray(
        [
            float(value)
            for line in lines[1:]
            for value in line.split()
        ],
        dtype=float,
    )

    nz = len(depths)

    if nz == 0 or values.size % nz != 0:
        raise ValueError(
            f"MOD数据量与深度层数不匹配："
            f"values={values.size}, nz={nz}"
        )

    nodes_per_layer = values.size // nz
    nxy = int(round(np.sqrt(nodes_per_layer)))

    if nxy * nxy != nodes_per_layer:
        raise ValueError(
            "当前转换脚本要求经纬度网格数量相同；"
            f"每层节点数={nodes_per_layer}，无法组成方形网格。"
        )

    model = values.reshape(nz, nxy, nxy)

    return depths, model


def read_inverted_model(path: Path):
    data = np.loadtxt(
        require_file(path, "T1最终反演模型"),
        dtype=float,
    )

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 4:
        raise ValueError(
            f"{path.name}必须至少有4列："
            "longitude latitude depth velocity"
        )

    data = data[:, :4]

    if not np.all(np.isfinite(data)):
        raise ValueError("T1反演模型中存在NaN或Inf。")

    keys = [
        (
            round(float(lon), 7),
            round(float(lat), 7),
            round(float(depth), 6),
        )
        for lon, lat, depth, _ in data
    ]

    if len(keys) != len(set(keys)):
        raise ValueError("T1反演模型中存在重复空间节点。")

    return data


def main():
    template_depths, template_model = read_mod_template(TEMPLATE_MOD)
    inverted = read_inverted_model(INPUT_MODEL)

    nz_template, ny_template, nx_template = template_model.shape

    lon_values = np.unique(inverted[:, 0])
    lat_values = np.unique(inverted[:, 1])
    depth_values = np.unique(inverted[:, 2])

    # DSurfTomo的i方向纬度由北向南递减；
    # j方向经度由西向东递增。
    lon_order = np.sort(lon_values)
    lat_order = np.sort(lat_values)[::-1]
    depth_order = np.sort(depth_values)

    nx_inner = len(lat_order)
    ny_inner = len(lon_order)
    nz_inner = len(depth_order)

    if nx_inner != nx_template - 2:
        raise ValueError(
            f"反演纬度节点数={nx_inner}，"
            f"但MOD内部应为{nx_template - 2}。"
        )

    if ny_inner != ny_template - 2:
        raise ValueError(
            f"反演经度节点数={ny_inner}，"
            f"但MOD内部应为{ny_template - 2}。"
        )

    output_model = template_model.copy()
    inserted_depths = []
    preserved_depths = []

    for depth in template_depths:
        if np.any(np.isclose(depth_order, depth, atol=1.0e-6)):
            inserted_depths.append(float(depth))
        else:
            preserved_depths.append(float(depth))

    for depth in depth_order:
        matching_template_depth = np.where(
            np.isclose(template_depths, depth, atol=1.0e-6)
        )[0]

        if len(matching_template_depth) != 1:
            raise ValueError(
                f"反演深度{depth:g} km在MOD模板中"
                f"匹配到{len(matching_template_depth)}层。"
            )

        k = int(matching_template_depth[0])
        interior = np.empty((ny_inner, nx_inner), dtype=float)

        for j, lon in enumerate(lon_order):
            for i, lat in enumerate(lat_order):
                row = inverted[
                    np.isclose(inverted[:, 0], lon, atol=1.0e-7)
                    & np.isclose(inverted[:, 1], lat, atol=1.0e-7)
                    & np.isclose(inverted[:, 2], depth, atol=1.0e-6)
                ]

                if len(row) != 1:
                    raise ValueError(
                        "缺少或重复节点："
                        f"lon={lon}, lat={lat}, depth={depth}"
                    )

                interior[j, i] = float(row[0, 3])

        output_model[k, 1:-1, 1:-1] = interior

    # 确认写入后的内部节点与原始反演文件逐点一致。
    maximum_insert_error = 0.0

    for depth in depth_order:
        k = int(
            np.where(
                np.isclose(template_depths, depth, atol=1.0e-6)
            )[0][0]
        )

        for j, lon in enumerate(lon_order):
            for i, lat in enumerate(lat_order):
                expected = inverted[
                    np.isclose(inverted[:, 0], lon, atol=1.0e-7)
                    & np.isclose(inverted[:, 1], lat, atol=1.0e-7)
                    & np.isclose(inverted[:, 2], depth, atol=1.0e-6),
                    3,
                ][0]

                actual = output_model[k, j + 1, i + 1]
                maximum_insert_error = max(
                    maximum_insert_error,
                    abs(float(actual) - float(expected)),
                )

    with OUTPUT_MOD.open("w", encoding="utf-8") as handle:
        handle.write(
            "  "
            + "  ".join(f"{depth:g}" for depth in template_depths)
            + "\n"
        )

        for k in range(nz_template):
            for j in range(ny_template):
                handle.write(
                    "  "
                    + "  ".join(
                        f"{velocity:8.5f}"
                        for velocity in output_model[k, j, :]
                    )
                    + "\n"
                )

    # 重新读取输出，检查数量和格式。
    output_lines = OUTPUT_MOD.read_text(
        encoding="utf-8",
        errors="ignore",
    ).splitlines()

    written_depths = np.asarray(
        [float(value) for value in output_lines[0].split()],
        dtype=float,
    )

    written_values = np.asarray(
        [
            float(value)
            for line in output_lines[1:]
            for value in line.split()
        ],
        dtype=float,
    )

    expected_value_count = nz_template * ny_template * nx_template

    if len(written_depths) != nz_template:
        raise RuntimeError("输出MOD的深度节点数量错误。")

    if written_values.size != expected_value_count:
        raise RuntimeError(
            f"输出MOD速度数量={written_values.size}，"
            f"应为{expected_value_count}。"
        )

    report = "\n".join([
        "DSurfTomo T1 → DSurfTomo-4D MOD 转换完成",
        "=" * 68,
        f"输入最终模型：{INPUT_MODEL.resolve()}",
        f"结构模板：{TEMPLATE_MOD.resolve()}",
        f"输出4D MOD：{OUTPUT_MOD.resolve()}",
        "",
        f"模板网格：{nx_template} × {ny_template} × {nz_template}",
        f"反演内部网格：{nx_inner} × {ny_inner} × {nz_inner}",
        f"输入反演节点数：{len(inverted)}",
        f"输出速度值数量：{written_values.size}",
        "",
        "写入的深度层："
        + ", ".join(f"{value:g}" for value in inserted_depths)
        + " km",
        "保留模板的深度层："
        + (
            ", ".join(f"{value:g}" for value in preserved_depths)
            if preserved_depths
            else "无"
        )
        + " km",
        "保留内容：最外一圈经纬度边界节点，以及未反演深度层。",
        f"写入节点最大核对误差：{maximum_insert_error:.3e} km/s",
        "",
        f"MOD模板 MD5：{md5(TEMPLATE_MOD)}",
        f"最终模型 MD5：{md5(INPUT_MODEL)}",
        f"输出MOD MD5：{md5(OUTPUT_MOD)}",
        "",
        "使用方法：",
        "将 MOD_4D 复制到4D实验目录并改名为 MOD。",
    ])

    REPORT_FILE.write_text(report + "\n", encoding="utf-8")

    print(report)


if __name__ == "__main__":
    main()
