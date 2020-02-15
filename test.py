#!/usr/bin/env python
# -*-coding: utf-8 -*-

"""
Usage: $ {1: program}.py
"""

import itertools

a = ["a", "b"]
b = ["c", "d", "e"]
c = [a]

print(list(itertools.product(*c)))
