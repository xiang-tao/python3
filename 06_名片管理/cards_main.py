#! /usr/bin/python3

import cards_tools
# 无限循环由用户决定什么时候退出，即输入0
while True:

    # 显示功能菜单
    cards_tools.show_menu()

    action_str = input("请选折希望执行的操作: ")
    print("您选择的操作是[%s]" % action_str)

    # 1 2 3 针对名片的操作
    if action_str in ["1", "2", "3"]:
        # 新增名片
        if action_str == "1":
            cards_tools.new_card()
        # 显示全部
        elif action_str == "2":
            cards_tools.show_all()
        # 查询名片
        elif action_str == "3":
            cards_tools.search_card()
        # 如果在程序开发时，不希望立刻编写分支内部的代码
        # 可以使用pass 关键字，表示一个占位符，能够保证程序的代码结构正确
        # 程序执行时，pass关键字不会执行任何操作
        # pass
        pass
    # 0退出系统
    elif action_str == "0":
        print("欢迎再次使用[名片管理系统]")
        break

    # 其他内容输入错误，需要提示用户
    else:
        print("您输入的不正确，请重新选择")
