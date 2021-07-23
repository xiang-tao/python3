def print_line(char, times):
    print(char * times)


def print_lines(char, times, cunt):
    """

    :param char:
    :param times:
    :param cunt:
    """
    for i in range(cunt):
        print_line(char, times)


print_lines("*", 30, 8)
