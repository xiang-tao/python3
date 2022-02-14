# 1.判断空白字符
# space_str = "  "
# space_str = "  \t\n\r"，空格，\t\n\r等字符认为是空白字符
space_str = "  a"
print(space_str.isspace())

# 2.判断字符串中是否只包含数字
num_str = "1"
# 1)都不能够判断小数
# num_str = "1.1"
# 2）unicode字符串
# num_str = "\u00b2"
# 3)中文字符
# num_str = "向涛"
print(num_str)
print(num_str.isdecimal())
print(num_str.isdigit())
print(num_str.isnumeric())
