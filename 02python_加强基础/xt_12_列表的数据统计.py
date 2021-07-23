name_list = ["张三", "李四", "王五", "张三"]

# len(leng 长度) 函数可以统计列表中元素的总数
list_len = len(name_list)
print("列表中一共有 %d 个元素" % list_len)

# count 方法可以统计列表中某一个数据出现的次数
count = name_list.count("张三")
print("张三出现了 %d 次" % count)
