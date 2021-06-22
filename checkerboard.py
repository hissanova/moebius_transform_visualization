from multiprocessing import Pool
import argparse
from itertools import product
from typing import Dict, List, NamedTuple, Tuple, Union
from dataclasses import dataclass

from numpy import exp, inf, pi, sqrt, tan

from circles import Num, Point2D, Line, CanonicalCircle, Circle, AppolonianCircle


class Range(NamedTuple):
    inf: Num
    sup: Num


@dataclass
class BoundedRegion:
    lower_bound: Circle
    upper_bound: Circle
    col_key: Num


def cot(tangent_angle: Num) -> Num:
    """It is assumed that foci are at (0,1) and (0,-1). TO-DO: generalize to arbitrary foci"""
    if tangent_angle < 0 or pi < tangent_angle:
        raise Exception("tangent angle a at focus must be between 0 <= a <= pi. Given value: {}".format(tangent_angle))
    if tangent_angle == 0:
        return + inf
    elif tangent_angle == pi:
        return - inf
    else:
        return 1/tan(tangent_angle)


def test_exceptions(values):
    for value in values:
        try:
            res = cot(value)
        except Exception as e:
            print(e)
        else:
            print(res)


test_exceptions([-1, 2*pi, 0, pi])
print(cot(0) > 0)
print(cot(pi) > 0)


def get_verti_bound(centre_x: Num) -> Union[Line, CanonicalCircle]:
    if centre_x is +inf:
        return Line(Point2D(-1, 0), 0)
    elif centre_x is -inf:
        return Line(Point2D(1, 0), 0)
    else:
        centre_x = float(centre_x)
        return CanonicalCircle(Point2D(centre_x, 0), sqrt(centre_x ** 2 + 1))


def mod_(x: Num, p: Num) -> Num:
    if 0 <= x and x <= p:
        return x
    elif x < 0:
        return mod_(x + p, p)
    else:
        return mod_(x - p, p)


def get_col_key(num: Num) -> Num:
    return num % 2


def get_vertical_bounded_regions(angle_list: List[Tuple[Num, Num]]) -> Tuple[List[Union[Line, CanonicalCircle]], List[BoundedRegion]]:
    verti_boundaries = []
    verti_region_grids = []
    for i, (lower_angle, upper_angle) in enumerate(angle_list):
        lower_centre_x = cot(lower_angle)
        lower_bound = get_verti_bound(lower_centre_x)
        verti_boundaries.append(lower_bound)
        upper_centre_x = cot(upper_angle)
        upper_bound = get_verti_bound(upper_centre_x)
        if isinstance(lower_bound, CanonicalCircle) and isinstance(upper_bound, CanonicalCircle):
            if upper_bound.c.x > lower_bound.c.x:
                upper_bound.flip_insideout()
        if isinstance(lower_bound, Line) and isinstance(upper_bound, CanonicalCircle):
            if lower_bound.v.x > 0:
                lower_bound.flip_insideout()
        verti_region_grids.append(BoundedRegion(lower_bound, upper_bound, get_col_key(i)))
    return verti_boundaries, verti_region_grids


def get_horizontal_bounded_regions(ratio_list: List[Tuple[Num, Num]],
                                   focus_list: List[Point2D]) -> Tuple[List[AppolonianCircle], List[BoundedRegion]]:
    focus1, focus2 = focus_list
    horiz_boundaries = []
    horiz_region_grids = []
    for i, (lower_ratio, upper_ratio) in enumerate(ratio_list):
        lower_bound = AppolonianCircle(focus1, focus2, lower_ratio)
        horiz_boundaries.append(lower_bound)
        upper_bound = AppolonianCircle(focus1, focus2, upper_ratio)
        horiz_region_grids.append(BoundedRegion(lower_bound, upper_bound, get_col_key(i)))
    return horiz_boundaries, horiz_region_grids


def get_angle_grids(angle_increment: Num = 0,
                    p: int = 6) -> List[Tuple[Num, Num]]:
    return [(mod_(i * pi/p + angle_increment, pi), mod_((i + 1) * pi/p + angle_increment, pi)) for i in range(p)]


def get_ratio_grids() -> List[Range]:
    val_list: List[Num]
    val_list = [exp(0.5 * i) for i in range(-4, 5)]
    val_list = [0.] + val_list + [inf]
    ratio_grids = []
    for i in range(len(val_list) - 1):
        ratio_grids.append(Range(val_list[i], val_list[i+1]))
    return ratio_grids


def make_checkerboard(angle_grids,
                      ratio_grids,
                      focus_list,
                      col_dict: Dict = {0: 'yellow', 1: 'lawngreen'}):
    # print(angle_grids)
    verti_boundaries, verti_region_grids = get_vertical_bounded_regions(angle_grids)
    # print(verti_region_grids)

    # print(ratio_grids)
    horiz_boundaries, horiz_region_grids = get_horizontal_bounded_regions(ratio_grids, focus_list)
    # print(horiz_region_grids)
    boundaries = verti_boundaries + horiz_boundaries

    region_list = []
    for horiz_region, verti_region in product(horiz_region_grids, verti_region_grids):
        if horiz_region.upper_bound.is_point:
            region0 = [horiz_region.lower_bound.alg_eq > 0]
        else:
            region0 = [horiz_region.lower_bound.alg_eq > 0, horiz_region.upper_bound.alg_eq < 0]
        # Plots two regions of symmetric difference of the two boundary circles: (not C_1 \cap C_2) \cup (C_1 and \cap C_2)
        region1 = [verti_region.lower_bound.alg_eq < 0, verti_region.upper_bound.alg_eq > 0]
        region2 = [verti_region.lower_bound.alg_eq > 0, verti_region.upper_bound.alg_eq < 0]
        col_key = get_col_key(horiz_region.col_key + verti_region.col_key)
        incol = col_dict[col_key]
        region_list.extend([(region1 + region0, incol), (region2 + region0, incol)])
    return region_list, boundaries


def render_objects(region_list,
                   boundaries,
                   x_range: Range,
                   y_range: Range,
                   point_list: List[Point2D] = [],
                   plot_points: int = 200,
                   bound_col: str = 'black'):
    # x_min, x_max = x_range
    # y_min, y_max = y_range
    # G = Graphics()
    # if len(point_list) > 0:
    #     for point in point_list:
    #         G += point2d(point)

    # for region, incol in region_list:
    #     G += region_plot(region, [x, x_min, x_max], [y, y_min, y_max], incol=incol, plot_points=plot_points)

    # for boundary in boundaries:
    #     region_eq = boundary.alg_eq
    #     G += implicit_plot(region_eq == 0 , [x, x_min, x_max], [y, y_min, y_max], color=bound_col, frame=False)
    # return G
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frame_num", type=int, default=40)
    parser.add_argument("--plot_points", type=int, default=200)
    return parser.parse_args()


args = parse_args()
y_range = Range(-6, 6)
x_range = Range(-6, 6)
focus1 = Point2D(0, 1)
focus2 = Point2D(0, -1)
focus_list = [focus1, focus2]
p = 6
plot_points = args.plot_points


def incremented_graph(increment: float):
    angle_grids = get_angle_grids(angle_increment=increment, p=p)
    ratio_grids = get_ratio_grids()
    region_list, boundaries = make_checkerboard(angle_grids, ratio_grids, focus_list)
    return render_objects(region_list,
                          boundaries,
                          x_range,
                          y_range,
                          point_list=focus_list,
                          plot_points=plot_points)


n = args.frame_num

param_list = [i/n * pi/3 for i in range(n)]

with Pool() as pool:
    frames = pool.map(incremented_graph, param_list)
# pool = ProcessPool()
# frames = list(map(incremented_graph, param_list))
# plt.show()

# a = animate(frames, xmin=x_range[0], xmax=x_range[1], ymin=y_range[0], ymax=y_range[1], figsize=[6, 6])
# a.save("moebius-transform-elliptic.gif")
