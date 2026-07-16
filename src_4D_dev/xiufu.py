fpath = "CalSurfG.f90"
with open(fpath, "r") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    # 为 Rayleigh 和 Love 波的所有赋值语句添加安全检查
    if "velf(1:nx*ny)=pvRc" in line:
        new_lines.append("          if (kmaxRc > 0) " + line.strip() + "\n")
    elif "velf(1:nx*ny)=pvRg" in line:
        new_lines.append("          if (kmaxRg > 0) " + line.strip() + "\n")
    elif "velf(1:nx*ny)=pvLc" in line:
        new_lines.append("          if (kmaxLc > 0) " + line.strip() + "\n")
    elif "velf(1:nx*ny)=pvLg" in line:
        new_lines.append("          if (kmaxLg > 0) " + line.strip() + "\n")
    else:
        new_lines.append(line)

with open(fpath, "w") as f:
    f.writelines(new_lines)
print("✅ CalSurfG.f90 安全逻辑加固完成！")
