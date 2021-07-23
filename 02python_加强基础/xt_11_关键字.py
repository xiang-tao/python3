name_list = ["张三", "李四", "王五"]

# 知道使用 del 关键字（delete）删除列表元素
# 提示：在日常开发中要从列表删除数据建议使用列表提供的方法
del name_list[1]

name = "小明"

del name
# 注意如果使用 del 关键词将变量从内存中删除
# 后续的代码就不能够使用这个变量了
# print(name)

print(name_list)
