def sum_1_num():
    """
    这是一个两个数字的求和函数
    :return:
    """
    num1 = 10
    num2 = 20
    result = num1 + num2
    print("%d + %d = %d" % (num1, num2, result))


sum_1_num()


def sum_2_num(a, b):
    """
    这是任意两个数字的求和函数
    :return:
    """
    result = a + b
    print("%d + %d = %d" % (a, b, result))


sum_2_num(1, 1)


def sum_3_num(a, b):
    """
    这是任意两个数字的求和函数，并使用return返回函数的数据
    return可以使得函数调用方也知道函数调用的结果（其实似乎没啥太大区别）
    其次return需要函数调用方使用变量来接收结果，不然不会将求和的结果显示
    并且return标志着函数的结束，如果接着在return的下方写代码那将变得毫无意义
    :return:
    """
    return a + b


sum_3_result = sum_3_num(1, 3)
print("计算结果： %d " % sum_3_result)
