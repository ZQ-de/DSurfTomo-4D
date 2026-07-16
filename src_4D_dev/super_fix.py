import re

with open('main.f90', 'r') as f:
    content = f.read()

# 修复那个被 & 折断的逻辑（确保不再多加 kmaxRc）
content = content.replace('kmaxRc+&', '&')
content = content.replace('kmaxRc+', '')

# 核心手术：找到分配内存的地方，给它增加 100 个格子的缓冲区
# 这样即便数据里有第 20、30 个周期的索引，也不会崩
content = content.replace('kmaxLc+kmaxRg+kmaxRc', '(kmaxLc+kmaxRg+kmaxRc+100)')

with open('main.f90', 'w') as f:
    f.write(content)

print("✅ 逻辑修正与内存缓冲区（+100）已设置完成！")
