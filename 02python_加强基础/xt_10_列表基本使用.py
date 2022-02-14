name_list = ["zhangsan", "lisi", "wangwu"]

# 1.取值和索引（index）
# IndexError: list index out of range - 列表索引超出范围
# print(name_list[3])
print(name_list[2])

# 知道数据内容查找数据位置
# ValueError: 'wangwuwu' is not in list - 该数据不在列表中，
# print(name_list.index("wangwuwu"))
print(name_list.index("wangwu"))

# 2.修改
# IndexError: list assignment index out of range - 列表指定的索引超出范围
# name_list[3] = "王小而"
name_list[1] = "李四"


# 3.增加
# append 方法可以向列表的末尾追加数据
name_list.append("王小二")

# insert 方法可以在列表指定的索引为位置插入数据
name_list.insert(1, "小美女")

# extend 方法可以把其他列表中的完整内容追加到当前列表的末尾
temp_list = ["sunwukong", "zhubajie", "shashidi"]
name_list.extend(temp_list)


# 4.删除
# remove 方法可以从列表中删除指定的数据
# 注意当有列表中有多个 wangwu 时候，remove 方法会删除列表中第一个 wangwu
name_list.remove("wangwu")

# pop 方法默认可以删除列表最后一个参数
name_list.pop()
# pop 方法可以指定要删除元素的索引，并删除该索引处的元素
name_list.pop(3)

# clear 方法可以清空列表
name_list.clear()


print(name_list)
