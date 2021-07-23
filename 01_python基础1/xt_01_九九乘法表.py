for i in range(9):
    for j in range(i+1):
        # 此处的end="",表示不换行，因为print()默认输出一次数据后换行，
        # \t是制表符（转义字符，比如\n表示换行符），作用就是使得垂直方向对齐。
        print("%d * %d = %d" % ((j + 1), (i + 1), ((j+1) * (i+1))), end="\t")
    print()
