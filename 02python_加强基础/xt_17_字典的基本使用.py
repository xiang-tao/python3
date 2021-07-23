# key:value
xiaoming_dict = {"name": "小明"}

# 1.取值 print(字典名[key])
print(xiaoming_dict["name"])
# 当name22(key)不再字典中时候会报错
# print(xiaoming_dict["name22"])

# 2.增加/修改
# 如果key不存在，会新增键值对
xiaoming_dict["age"] = 18
# 如果key存在，会修改已存在的键值对
xiaoming_dict["name"] = "小红"

# 3.删除
# 在字典的删除方法中没有提供remove，而是pop方法，
# 即 字典名.pop("key")即可删除对应的value
xiaoming_dict.pop("name")
# 当key不存在时就会报错
# xiaoming_dict.pop("name11")

print(xiaoming_dict)
