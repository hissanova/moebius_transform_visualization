from typing import Callable, List, NamedTuple, Union
from abc import ABC, abstractmethod
from numpy import inf
from sympy import var

Num = Union[float, int]

x, y = var("x y")


class Point2D(NamedTuple):
    x: Num
    y: Num


class CharacteristicFunction(ABC):
    def __init__(self,
                 f: Callable[[Point2D], Num]):
        self.f = f

    def __call__(self, point: Point2D) -> int:
        return int(self.f(point) > 0)

    def __add__(self, other: 'CharacteristicFunction') -> 'CharacteristicFunction':
        def _func(p: Point2D) -> int:
            return (self(p) + other(p)) % 2
        return _func


class Line:
    def __init__(self,
                 orthonormal: Point2D,
                 d: Num):
        self.v = orthonormal
        self.d = d
        self._insideout = False

    def flip_insideout(self) -> None:
        self._insideout = not self._insideout

    def get_alg_eq(self):
        if self._insideout:
            return -(self.v.x * x + self.v.y * y - self.d)
        else:
            return self.v.x * x + self.v.y * y - self.d


class Circle:
    def __init__(self,
                 centre: Point2D,
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
            return -((x - self.c.x) ** 2 + (y - self.c.x) ** 2 - self.r ** 2)
        else:
            return (x - self.c.x) ** 2 + (y - self.c.y) ** 2 - self.r ** 2


class AppolonianCircle:
    def __init__(self,
                 focus1: Point2D,
                 focus2: Point2D,
                 ratio: Num):
        self.f1 = focus1
        self.f2 = focus2
        self.r = ratio

    def is_point(self):
        return self.r == 0 or self.r is inf

    def get_alg_eq(self):
        if self.r == 0:
            return (x - self.f1.x) ** 2 + (y - self.f1.y) ** 2
        elif self.r is inf:
            return (x - self.f2.x) ** 2 + (y - self.f2.y) ** 2
        else:
            return (x - self.f1.x) ** 2 + (y - self.f1.y) ** 2 - self.r * ((x - self.f2.x) ** 2 + (y - self.f2.y) ** 2)
