xiaoming_dict = {"name": "小明",
                 "age": 18}

# 1.统计键值对数量
print(len(xiaoming_dict))

# 2.合并字典
temp_dict = {"height": 1.75,
             "age": 20}
# 注意在update方法中如果包含已存在的键值对则会覆盖，
# 例如这里的age，没有则会新增，例如height
xiaoming_dict.update(temp_dict)

# 3.清空字典
xiaoming_dict.clear()

print(xiaoming_dict)
