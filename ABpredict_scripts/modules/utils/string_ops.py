#! /usr/bin/env python3
"""
commonly used list related operations
"""
import random
import string


def rand_str(size: int) -> str:
    """rand_str
    returns a random string of digits and characters, length size
    :param size: size of string to produce
    :type size: int
    :rtype: str
    """
    return ''.join(random.choices(string.ascii_lowercase +
                                  string.digits, k=size))
