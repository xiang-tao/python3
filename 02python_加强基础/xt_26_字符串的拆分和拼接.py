# 假设以下内容来源于网络抓取
# 要求：
# 1.将字符串中的空白字符全部去掉
# 2.再使用" "作为分隔符，拼接成一个整齐的字符串
poem_str = "\t\n登鹤雀楼  王之涣  \r白日依山尽 \t 黄河入海流\n\r 欲穷千里目  \n更上一层楼"


print(poem_str)
# 1.拆分字符串
poem_list = poem_str.split()
print(poem_list)

# 2.合并字符串
# result = "--".join(poem_list)
result = "  ".join(poem_list)
print(result)
