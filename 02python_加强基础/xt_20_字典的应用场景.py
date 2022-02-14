# 使用 多个键值对存储描述一个物体的相关信息 - 描述更复杂
# 的数据信息，将多个字典放在一个列表中再进行遍历

card_list = [
    {"name": "xiaoming",
     "qq": 123456,
     "phone": 10086},
    {"name": "lisi",
     "qq": 123456,
     "phone": 1008611}
]

for card_info in card_list:
    print(card_info)
