name = "小明"
# python 解释器知道下方定义了一个函数


def say_hello():
    """
    这是对函数作用说明注释的地方，例如当下面的代码调用了该函数想要
    查看这些注释时候可以点击函数名后按下 crtl+Q即可
    :return: 填写返回值
    """
    print('hello python')
    print('hello matlab')
    print('hello world')


print(name)
# 只有在程序中主动调用函数才会让函数执行
say_hello()

print(name)
