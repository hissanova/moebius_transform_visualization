from typing import List, Tuple, Union

from sympy import var

Num = Union[float, int]

x, y = var("x y")


class Line:
    def __init__(self,
                 orthonormal: List[float],
                 d: Num):
        self.v = orthonormal
        self.d = d
        self._insideout = False

    def flip_insideout(self) -> None:
        self._insideout = not self._insideout

    def get_alg_eq(self):
        if self._insideout:
            return -(self.v[0] * x + self.v[1] * y - self.d)
        else:
            return self.v[0] * x + self.v[1] * y - self.d


class Circle:
    def __init__(self,
                 centre: Num,
                 radius: Num,
                 insideout: bool = False):
        self.c = centre
        self.r = radius
        self._insideout = insideout

    def is_isideout(self) -> bool:
        return self._insideout

    def flip_insideout(self) -> None:
        self._insideout = not self._insideout

    def is_point(self) -> bool:
        return self.r == 0

    def get_alg_eq(self):
        if self._insideout:
            return -((x - self.c[0]) ^ 2 + (y - self.c[1]) ^ 2 - self.r ^ 2)
        else:
            return (x - self.c[0]) ^ 2 + (y - self.c[1]) ^ 2 - self.r ^ 2


class AppolonianCircle:
    def __init__(self, focus1, focus2, ratio):
        self.f1 = focus1
        self.f2 = focus2
        self.r = ratio

    def is_point(self):
        return self.r == 0 or self.r is Infinity

    def get_alg_eq(self):
        if self.r == 0:
            return (x - self.f1[0]) ^ 2 + (y - self.f1[1]) ^ 2
        elif self.r is Infinity:
            return (x - self.f2[0]) ^ 2 + (y - self.f2[1]) ^ 2
        else:
            return (x - self.f1[0]) ^ 2 + (y - self.f1[1]) ^ 2 - self.r * ((x - self.f2[0]) ^ 2 + (y - self.f2[1]) ^ 2)
