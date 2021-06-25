from numpy import inf, pi, tan

from src.circles import Num


def cot(tangent_angle: Num) -> Num:
    """It is assumed that foci are at (0,1) and (0,-1). TO-DO: generalize to arbitrary foci"""
    if tangent_angle < 0 or pi < tangent_angle:
        raise Exception(
            f"tangent angle a at focus must be between 0 <= a <= pi. Given value: {tangent_angle}"
        )
    if tangent_angle == 0:
        return +inf
    elif tangent_angle == pi:
        return -inf
    else:
        return 1 / tan(tangent_angle)


def mod_(x: Num, p: Num) -> Num:
    if 0 <= x and x <= p:
        return x
    elif x < 0:
        return mod_(x + p, p)
    else:
        return mod_(x - p, p)
