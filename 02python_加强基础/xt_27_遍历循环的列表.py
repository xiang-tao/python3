students = [{"name": "向涛"},
            {"name": "王唯"},
            {"name": "赵世杰"},
            {"name": "赵浩科"}]

find_name = "老王"
for stu_dict in students:
    print(stu_dict)

    if stu_dict["name"] == find_name:
        print("找到了　%s" % find_name)
        break
else:
    print("抱歉没有找到最 %s " % find_name)
print("循环结束")
