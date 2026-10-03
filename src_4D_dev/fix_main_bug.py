import re

with open('main.f90', 'r') as f:
    content = f.read()

# 这种匹配方式可以无视换行符和 & 符号，精准定位那行错误逻辑
pattern = r"if\(wavetp\.eq\.1\.and\.veltp\.eq\.1\)\s+knum=kmaxLc\+kmaxRg\+&\s+kmaxRc\+period"
replacement = "if(wavetp.eq.1.and.veltp.eq.1) knum=kmaxLc+kmaxRg+period"

if re.search(pattern, content):
    new_content = re.sub(pattern, replacement, content)
    with open('main.f90', 'w') as f:
        f.write(new_content)
    print("成功：已修复 main.f90 中的换行 Bug！")
else:
    # 如果没找到换行版的，再尝试寻找单行版的
    pattern2 = r"knum=kmaxLc\+kmaxRg\+kmaxRc\+period"
    replacement2 = "knum=kmaxLc+kmaxRg+period"
    if re.search(pattern2, content):
        new_content = re.sub(pattern2, replacement2, content)
        with open('main.f90', 'w') as f:
            f.write(new_content)
        print("成功：已修复 main.f90 中的单行 Bug！")
    else:
        print("警告：未找到匹配的错误代码行，请检查文件内容。")
