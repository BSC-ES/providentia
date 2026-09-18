""" Functions for plot formatting """

import copy
from calendar import monthrange
from datetime import datetime, timedelta
import os

import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.dates import date2num, num2date
from matplotlib.figure import Figure
from matplotlib import ticker
from matplotlib import transforms as mtransforms
import numpy as np
import pandas as pd
from PIL import Image

from providentia.auxiliar import CURRENT_PATH, join, get_map_colours
from .plot_aux import get_land_polygon_resolution, set_map_extent
from .read_aux import get_periodic_relevant_temporal_resolutions
from .plot_options import (
    annotation,
    model_domain,
    linear_regression,
    log_axes,
    smooth,
    threshold,
)
from .statistics import get_z_statistic_info
from .warnings_prv import show_message

Image.MAX_IMAGE_PIXELS = None

# matplotlib's default font (DejaVu Sans) has no CJK glyphs, so Chinese
# station names/metadata render as boxes. Add whichever CJK font is already
# on the machine as a fallback, keeping DejaVu Sans as the primary font.
# Fixed install paths are probed as well as matplotlib's font cache, as the
# frozen Mac app (bin/Mac/Providentia) does not find the OS fonts on its own.
_CJK_FONT_CANDIDATES = [
    # macOS
    "PingFang SC", "PingFang TC", "PingFang HK", "Hiragino Sans GB",
    "Heiti SC", "Heiti TC", "STHeiti", "Songti SC", "Apple SD Gothic Neo",
    # Linux
    "Noto Sans CJK SC", "Noto Sans CJK TC", "Noto Sans CJK JP",
    "Noto Sans SC", "Noto Sans TC", "WenQuanYi Zen Hei",
    "WenQuanYi Micro Hei", "Droid Sans Fallback", "Source Han Sans SC",
    "Source Han Sans TC",
    # Windows
    "Microsoft YaHei", "Microsoft JhengHei", "SimHei", "SimSun",
]

# fixed, well-known install paths for the same fonts, per OS - tried
# directly against the filesystem, bypassing font_manager's own
# (potentially incomplete, e.g. inside a frozen app bundle) directory scan
_CJK_FONT_PATHS = [
    # macOS
    "/System/Library/Fonts/PingFang.ttc",
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/STHeiti Light.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/Supplemental/Songti.ttc",
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
    "/Library/Fonts/Arial Unicode.ttf",
    # Linux
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJKsc-Regular.otf",
    "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
    "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
    "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf",
    "/usr/share/fonts/truetype/arphic/uming.ttc",
    # Windows
    r"C:\Windows\Fonts\msyh.ttc",
    r"C:\Windows\Fonts\simhei.ttf",
    r"C:\Windows\Fonts\simsun.ttc",
    r"C:\Windows\Fonts\msjh.ttc",
]


def enable_cjk_font_fallback():
    """
    Make any already-installed CJK-capable font(s) available as a
    fallback for CJK text, without changing which font non-CJK (the
    vast majority of the dashboard's) text uses.

    This sets rcParams['font.family'] itself, deliberately not just
    rcParams['font.sans-serif']: matplotlib only does real per-glyph
    fallback across *multiple* concrete fonts when font.family holds a
    list of concrete font names. Left at its default generic alias
    value (['sans-serif']) and with only font.sans-serif edited,
    findfont() just resolves that alias to a single, primary concrete
    font - the first available name in font.sans-serif - and uses that
    one font for everything. An earlier version of this fix prepended
    the CJK font there, which "worked" for CJK glyphs but also made
    that CJK font (which incidentally covers Latin/ASCII fine too) the
    primary font for *all* text, visibly replacing DejaVu Sans
    dashboard-wide - a real, user-visible regression, not the intended
    surgical fallback. Setting font.family directly avoids that: the
    existing primary font (DejaVu Sans, matplotlib's bundled default)
    stays first and so stays primary for everything it already covers;
    the CJK candidates are only ever reached for glyphs it's missing.

    Not just called once at import time here: plot_aux.py calls
    seaborn's sns.reset_orig() on every plot-parameters refresh, which
    restores *all* rcParams (mpl.rcParams.update(mpl.rcParamsOrig)) to
    their state from before this module ever ran, silently wiping this
    fallback back out mid-session. So Plotting.make_metadata() also
    calls this again right before it draws text, to reinstate it every
    time regardless of what reset happened in between. It's cheap
    (a handful of dict lookups plus a handful of os.path.isfile() checks)
    so re-running it per metadata draw is not a concern.
    """
    found = []

    # first: whatever the environment's own font cache already knows about
    available = {f.name for f in font_manager.fontManager.ttflist}
    found += [name for name in _CJK_FONT_CANDIDATES if name in available]

    # second: explicit, fixed OS install paths, in case the font cache
    # missed them (e.g. a frozen app bundle's incomplete font scan) -
    # registers the file directly with font_manager rather than relying
    # on it having found the font on its own
    for path in _CJK_FONT_PATHS:
        if os.path.isfile(path):
            try:
                font_manager.fontManager.addfont(path)
                name = font_manager.FontProperties(fname=path).get_name()
            except Exception:
                continue
            if name not in found:
                found.append(name)

    if found:
        # font.sans-serif itself is left untouched throughout, so its
        # first entry is always the existing, unmodified primary font
        # (DejaVu Sans by default) - used here as font.family's own
        # primary, with the CJK fonts appended purely as fallback
        primary = mpl.rcParams["font.sans-serif"][0]
        mpl.rcParams["font.family"] = [primary] + [
            name for name in found if name != primary
        ]


enable_cjk_font_fallback()


def set_equal_axes(ax, plot_options, plot_characteristics, base_plot_type):
    """
    Set equal aspect ratio and axis limits for a plot (useful for scatter plots).

    Parameters
    ----------
    ax : object
        Matplotlib axis to set equal axes.
    plot_options : list
        Active plot options.
    plot_characteristics : dict
        Plot characteristics including optional 'xlim' and 'ylim'.
    base_plot_type : str
        Base plot type (without statistical overlays).
    """

    # if only one axis is on a log scale, then set box aspect to 1 to get square shape
    if (("logx" in plot_options) and ("logy" not in plot_options)) or (
        ("logx" not in plot_options) and ("logy" in plot_options)
    ):
        ax.set_box_aspect(1)
    else:
        # set equal aspect
        ax.set_aspect(aspect="equal", adjustable="box")

    if len(ax.lines) == 0 and base_plot_type == "scatter":
        return None

    if ("xlim" not in plot_characteristics) and ("ylim" not in plot_characteristics):
        # get min and max values for axes from plotted data
        xmin, xmax = get_data_lims(ax, "xlim", plot_options)
        ymin, ymax = get_data_lims(ax, "ylim", plot_options)

        # compare min and max lims across axes
        if xmin < ymin:
            axmin = xmin
        else:
            axmin = ymin
        if xmax > ymax:
            axmax = xmax
        else:
            axmax = ymax

        # set equal lims
        ax.set_xlim(axmin, axmax)
        ax.set_ylim(axmin, axmax)
    elif ("xlim" not in plot_characteristics) and ("ylim" in plot_characteristics):
        # set xlim as ylim if ylim is passed
        if isinstance(plot_characteristics["ylim"], dict):
            ax.set_xlim(**plot_characteristics["ylim"])
        else:
            ax.set_xlim(plot_characteristics["ylim"])
    elif ("xlim" in plot_characteristics) and ("ylim" not in plot_characteristics):
        # set ylim as ylim if xlim is passed
        if isinstance(plot_characteristics["xlim"], dict):
            ax.set_ylim(**plot_characteristics["xlim"])
        else:
            ax.set_ylim(plot_characteristics["xlim"])
    elif ("xlim" in plot_characteristics) and ("ylim" in plot_characteristics):
        # set both limits
        if isinstance(plot_characteristics["xlim"], dict):
            ax.set_xlim(**plot_characteristics["xlim"])
        else:
            ax.set_xlim(plot_characteristics["xlim"])
        if isinstance(plot_characteristics["ylim"], dict):
            ax.set_ylim(**plot_characteristics["ylim"])
        else:
            ax.set_ylim(plot_characteristics["ylim"])


# fixed hierarchy of calendar-aligned tick steps, finest to coarsest (hour
# multiples are every divisor of 24, so every sub-daily grid puts a tick on
# midnight; day multiples are fine-grained to give the search a close density
# match, and are phased on a month start so the 1st carries a label)
# how many passes the view is allowed to grow by to fit its end labels, and
# the overhang in pixels below which it is left alone
_EDGE_LABEL_FIT_PASSES = 5
_EDGE_LABEL_FIT_TOLERANCE = 0.5

# months to look through for a month start inside a view
_MONTH_START_SEARCH_LIMIT = 24

# tick kinds landing on dates a reader recognises, rather than on a stride
# counted from wherever the data happens to start
_CALENDAR_TICK_KINDS = (
    "fifthmonth",
    "thirdmonth",
    "semimonth",
    "month",
    "year",
)

# the days of the month each grid lands on, coarsest first - the order they
# are tried in, so the evenest grid that can fill the view is the one used
_MONTH_DAY_GRIDS = {
    "semimonth": (1, 15),
    "thirdmonth": (1, 10, 20),
    "fifthmonth": (1, 5, 10, 15, 20, 25),
}

# a grid has to put at least this many labels in the view to be worth using,
# below which the next finer one is tried and, failing all of them, a plain
# stride counted from the data
_MIN_MONTH_GRID_TICKS = 3

_TIMESERIES_TICK_STEPS = [
    ("hour", 1),
    ("hour", 2),
    ("hour", 3),
    ("hour", 4),
    ("hour", 6),
    ("hour", 8),
    ("hour", 12),
    ("hour", 24),
    # no day stride longer than a fortnight: past that the month and year
    # tiers cover the same spans one label per month, where a 45- or 90-day
    # stride lands on dates of no significance
    ("day", 1), ("day", 2), ("day", 3), ("day", 4), ("day", 5),
    ("day", 6), ("day", 7), ("day", 8), ("day", 9), ("day", 10),
    ("day", 12), ("day", 14),
    ("fifthmonth", 1),
    ("thirdmonth", 1),
    ("semimonth", 1),
    ("month", 1),
    ("month", 3),
    ("month", 6),
    ("year", 1),
    ("year", 2),
    ("year", 5),
    ("year", 10),
    ("year", 25),
    ("year", 50),
    ("year", 100),
    ("year", 250),
    ("year", 500),
    ("year", 1000),
]


def _parse_yyyymmdd(value):
    """
    Parse a configured start/end date (an int or numeric string like 20180101)
    as a datetime. Callers treat None as "the true loaded range isn't known
    here", not an error.

    Parameters
    ----------
    value : int or str or None
        Date to parse

    Returns
    -------
    datetime.datetime or None
        Parsed date, or None if missing/unparseable
    """
    if value is None:
        return None
    try:
        text = str(int(value))
        return datetime(int(text[:4]), int(text[4:6]), int(text[6:8]))
    except (TypeError, ValueError):
        return None


# named temporal resolutions in seconds, used to recognise a view zoomed down
# to a single observation, so it can be labelled with that observation's own
# time rather than the two edges of the span
_RESOLUTION_SECONDS = {
    "hourly": 3600,
    "hourly_instantaneous": 3600,
    "3hourly": 3 * 3600,
    "3hourly_instantaneous": 3 * 3600,
    "6hourly": 6 * 3600,
    "6hourly_instantaneous": 6 * 3600,
    "daily": 86400,
}


def _format_single_observation_tick(dt):
    """
    Get label text for the one tick shown when a view has been zoomed down to
    at most a single observation. Full year-month-day precision always, with
    hour/minute/second appended only as far as the timestamp needs.

    Parameters
    ----------
    dt : datetime.datetime
        Time of the observation

    Returns
    -------
    str
        Label text
    """
    if (dt.hour, dt.minute, dt.second, dt.microsecond) == (0, 0, 0, 0):
        return dt.strftime("%Y-%m-%d")
    if (dt.minute, dt.second, dt.microsecond) == (0, 0, 0):
        return dt.strftime("%Y-%m-%d %Hh")
    if (dt.second, dt.microsecond) == (0, 0):
        return dt.strftime("%Y-%m-%d %H:%M")
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def _month_start_within(left, right):
    """
    The first month start inside the view, if it holds one.

    What a day stride counts from, so that the 1st of the month carries a
    label rather than the grid running through it at whatever offset the
    data's own first day happens to give.

    Parameters
    ----------
    left : datetime.datetime
        Start of the range
    right : datetime.datetime
        End of the range

    Returns
    -------
    datetime.datetime or None
        The first month start in the range, or None if it holds none
    """

    year, month = left.year, left.month
    for _ in range(_MONTH_START_SEARCH_LIMIT):
        current = datetime(year, month, 1)
        if current > right:
            return None
        if current >= left:
            return current
        year, month = (year + 1, 1) if month == 12 else (year, month + 1)

    return None


def _aligned_timeseries_ticks(left, right, kind, multiple, anchor=None):
    """
    Get all datetimes in [left, right] landing exactly on a calendar boundary
    for the given step, e.g. kind="hour", multiple=3 gives every 3rd hour on
    the clock (00:00, 03:00, ...), never an arbitrary offset grid.

    Parameters
    ----------
    left : datetime.datetime
        Start of the range
    right : datetime.datetime
        End of the range
    kind : str
        Step kind ("hour", "day", "semimonth", "month" or "year")
    multiple : int
        Step multiple
    anchor : datetime.datetime, optional
        Day the multi-day grid counts from, normally the first day of the
        loaded data - so that the data's own start and end land on the grid
        and can be labelled, and so the ticks stay put as the view is panned.
        Left out, the grid counts from an absolute day ordinal instead.

    Returns
    -------
    list
        Aligned tick datetimes
    """

    ticks = []

    if kind == "hour":
        current = left.replace(minute=0, second=0, microsecond=0)
        if current < left:
            current += timedelta(hours=1)
        remainder = current.hour % multiple
        if remainder:
            current += timedelta(hours=multiple - remainder)
        step = timedelta(hours=multiple)
        while current <= right:
            ticks.append(current)
            current += step

    elif kind == "day":
        # counted from the data's own first day where one is known, so that
        # the start and end of the loaded range fall on the grid and get
        # labelled; otherwise day-ordinal aligned (day 1 is 0001-01-01), so a
        # multi-day stride still lands on a fixed set of calendar days.
        # Only ever used on views too short for the month grid, which is why
        # it can count from the data rather than from a month: over a few
        # days, dates a fixed stride apart are easy enough to follow
        start_date = left.date()
        if (left.hour, left.minute, left.second, left.microsecond) != (
            0,
            0,
            0,
            0,
        ):
            start_date += timedelta(days=1)
        ordinal = start_date.toordinal()
        origin = anchor.date().toordinal() if anchor is not None else 0
        remainder = (ordinal - origin) % multiple
        if remainder:
            ordinal += multiple - remainder
        step = timedelta(days=multiple)
        current = datetime.fromordinal(ordinal)
        while current <= right:
            ticks.append(current)
            current += step

    elif kind in _MONTH_DAY_GRIDS:
        # fixed days of every month in range - the 1st and the 15th, or the
        # 1st, 10th and 20th, or every fifth day. Dates a reader tracks by,
        # at the cost of one gap per month running long, as a month divides
        # into equal parts only by luck: the coarser the grid the smaller
        # that cost, which is why the coarsest that fills the view is the one
        # used. `multiple` is unused, but kept for a consistent per-kind
        # signature
        year, month = left.year, left.month
        while True:
            for day in _MONTH_DAY_GRIDS[kind]:
                if day > monthrange(year, month)[1]:
                    continue
                current = datetime(year, month, day)
                if current > right:
                    ticks.sort()
                    return ticks
                if current >= left:
                    ticks.append(current)
            year, month = (year + 1, 1) if month == 12 else (year, month + 1)

    elif kind == "month":
        # months counted as a single index from January of year 0, so a
        # step of e.g. 3 always lands on Jan/Apr/Jul/Oct, never an
        # offset depending on where the visible range happens to start
        month_index = left.year * 12 + (left.month - 1)
        if (left.day, left.hour, left.minute, left.second, left.microsecond) != (
            1,
            0,
            0,
            0,
            0,
        ):
            month_index += 1
        remainder = month_index % multiple
        if remainder:
            month_index += multiple - remainder
        while True:
            year, month = divmod(month_index, 12)
            current = datetime(year, month + 1, 1)
            if current > right:
                break
            ticks.append(current)
            month_index += multiple

    elif kind == "year":
        year = left.year
        if (left.month, left.day, left.hour, left.minute, left.second) != (
            1,
            1,
            0,
            0,
            0,
        ):
            year += 1
        remainder = year % multiple
        if remainder:
            year += multiple - remainder
        current = datetime(year, 1, 1)
        while current <= right:
            ticks.append(current)
            current = current.replace(year=current.year + multiple)

    return ticks


def _format_timeseries_tick(dt, kind):
    """
    Get label text for one aligned tick, precise enough for its own step and
    no more. "hour" ticks always carry their date too, so each reads
    unambiguously without depending on a nearby label for its date.

    Parameters
    ----------
    dt : datetime.datetime
        Tick time
    kind : str
        Step kind of the tier the tick belongs to

    Returns
    -------
    str
        Label text
    """

    if kind == "year":
        return dt.strftime("%Y")
    if kind == "month":
        return dt.strftime("%Y-%m")
    if kind in ("day",) or kind in _MONTH_DAY_GRIDS:
        return dt.strftime("%Y-%m-%d")
    return dt.strftime("%Y-%m-%d %Hh")


def _plotted_time_extent(ax):
    """
    Get the first and last time of the data plotted on a timeseries axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Timeseries axis

    Returns
    -------
    tuple of datetime.datetime or None
        (first, last) plotted time, or None if nothing is plotted in data
        coordinates
    """

    lines = [
        line
        for line in ax.lines
        if line.get_transform().contains_branch_seperately(ax.transData)[0]
        and len(line.get_xdata())
    ]
    visible_lines = [line for line in lines if line.get_visible()]
    lines = visible_lines or lines
    if not lines:
        return None

    x_mins, x_maxs = [], []
    for line in lines:
        x = np.asarray(ax.convert_xunits(line.get_xdata()), dtype=float)
        if np.all(np.isnan(x)):
            continue
        x_mins.append(np.nanmin(x))
        x_maxs.append(np.nanmax(x))
    if not x_mins:
        return None

    return (
        mpl.dates.num2date(min(x_mins)).replace(tzinfo=None),
        mpl.dates.num2date(max(x_maxs)).replace(tzinfo=None),
    )


def _format_free_timeseries_ticks(dates):
    """
    Get label texts for ticks not aligned to a calendar step (manual "n_ticks"
    or forced edges, both on whole hours), all at the precision the least
    round of them needs.

    Parameters
    ----------
    dates : list of datetime.datetime
        Tick times

    Returns
    -------
    list of str
        Label texts
    """

    if all(d.hour == 0 for d in dates):
        fmt = "%Y-%m-%d"
    else:
        fmt = "%Y-%m-%d %Hh"
    return [d.strftime(fmt) for d in dates]


def _set_timeseries_tick_alignment(tick_labels):
    """
    Centre every tick label on its tick, and reset each label's transform back
    to its un-nudged state. Both are set explicitly even when already correct,
    as matplotlib reuses a pool of Text objects across set_ticks() calls, so a
    label can otherwise inherit an alignment or a nudge left on a reused object
    by a previous call and quietly bias this call's own measurements.

    Parameters
    ----------
    tick_labels : list
        Tick label Text objects
    """

    for label in tick_labels:
        label.set_ha("center")
        base_transform = getattr(label, "_ptv_base_transform", None)
        if base_transform is None:
            base_transform = label.get_transform()
            label._ptv_base_transform = base_transform
        label.set_transform(base_transform)


def _edge_tick_text(dt, kind):
    """
    Get one edge's label text at the chosen tier - normally the same format as
    any interior tick of that tier, so the edge reads like it belongs on the
    same axis. "month", "year" and "hour" formats drop precision that is right
    for an edge genuinely on that tier's boundary but would hide a real
    difference otherwise, so each falls back one precision level finer for an
    edge that is not on its own boundary.

    Parameters
    ----------
    dt : datetime.datetime
        Edge time
    kind : str
        Step kind of the chosen tier

    Returns
    -------
    str
        Label text
    """
    if kind == "year" and (dt.month, dt.day) != (1, 1):
        return dt.strftime("%Y-%m-%d")
    if kind == "month" and dt.day != 1:
        return dt.strftime("%Y-%m-%d")
    if kind == "hour" and (dt.minute, dt.second, dt.microsecond) != (0, 0, 0):
        if dt.second or dt.microsecond:
            return dt.strftime("%Y-%m-%d %H:%M:%S")
        return dt.strftime("%Y-%m-%d %H:%M")
    return _format_timeseries_tick(dt, kind)


def _timeseries_labels_fit(ax, renderer, dates, texts, min_gap_pixels, xlim):
    """
    Whether every one of these labels can be shown at once without crowding
    its neighbours.

    Installed as real ticks and measured as rendered, rather than predicted
    from font metrics: a predicted width is only as good as the font and DPI
    assumptions behind it, and those differ between wherever this is tested
    and wherever it runs (font substitution, HiDPI scaling, hinting). The
    renderer is reused as-is rather than redrawing the figure per candidate,
    as get_window_extent() lays each string out against it afresh.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to install the ticks on
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure the labels against
    dates : list
        Candidate tick datetimes
    texts : list
        Their label text
    min_gap_pixels : float
        Smallest gap allowed between two neighbouring labels
    xlim : tuple
        Numeric x-axis limits to restore after setting ticks

    Returns
    -------
    bool
        Whether they all fit
    """

    ax.xaxis.set_ticks(dates, labels=texts)
    # matplotlib silently expands the axis's data limits to fit any tick
    # outside the current view, even with autoscale off, which would corrupt
    # every later measurement
    ax.set_xlim(*xlim)
    tick_labels = ax.xaxis.get_majorticklabels()
    _set_timeseries_tick_alignment(tick_labels)

    boxes = []
    for label in tick_labels:
        bbox = label.get_window_extent(renderer)
        # a small safety margin around the measured box, rather than trusting
        # it to the last pixel - two labels sitting exactly min_gap_pixels
        # apart with nothing to spare would still read as touching
        pad = (bbox.width * 0.15) / 2
        boxes.append((bbox.x0 - pad, bbox.x1 + pad))

    boxes.sort()
    for (_, previous_right), (next_left, _) in zip(boxes, boxes[1:]):
        if next_left - previous_right < min_gap_pixels:
            return False

    return True


def _widen_for_edge_labels(ax, renderer, dates):
    """
    Give the axis enough room either side for its outermost labels to sit
    centred under their ticks.

    A label centred on a tick at the very edge of the view has half of itself
    hanging outside the axis, and was previously shunted sideways to bring it
    back in - which reads as an off-centre label at one end of an otherwise
    even row. It is left centred instead, and the view stretched only if it
    would otherwise run off the figure entirely.

    Hanging over the edge of the axis is measured against the figure rather
    than the axis itself: a tick label is drawn outside its axis in any case,
    and demanding it fit within the axis box widened a narrow panel's view by
    a quarter to make room for a date - which then left too little room
    between the ticks for the labels that had just been chosen, collapsing a
    forty-day view to two labels.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis holding the ticks
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure the labels against
    dates : list
        The tick datetimes that were installed
    """

    if not dates:
        return None

    tick_labels = ax.xaxis.get_majorticklabels()
    if not tick_labels:
        return None

    axis_box = ax.get_window_extent(renderer)
    if axis_box.width <= 0:
        return None

    # what the labels must stay inside
    container = ax.figure.bbox

    # widening the view spreads the same pixels over a longer span, which
    # pulls the end ticks inwards - so the room needed is a little less than
    # the overhang measured before it. Rather than solve that, the step is
    # repeated until nothing reaches past the axis any more, which takes two
    # or three passes; the shrinking correction each time makes it converge
    for _ in range(_EDGE_LABEL_FIT_PASSES):
        left_limit, right_limit = ax.get_xlim()
        span = right_limit - left_limit
        if span <= 0:
            return None

        units_per_pixel = span / axis_box.width
        first_box = tick_labels[0].get_window_extent(renderer)
        last_box = tick_labels[-1].get_window_extent(renderer)
        left_overhang = max(0.0, container.x0 - first_box.x0)
        right_overhang = max(0.0, last_box.x1 - container.x1)

        if max(left_overhang, right_overhang) < _EDGE_LABEL_FIT_TOLERANCE:
            break

        ax.set_xlim(
            left_limit - (left_overhang * units_per_pixel),
            right_limit + (right_overhang * units_per_pixel),
        )

    return None


def _centre_view_on(ax, moment, span_seconds):
    """
    Put one instant in the middle of the axis, with the same width of view
    either side of it.

    Used where there is a single observation to show: left where it fell, it
    can sit anywhere across the panel - hard against one edge if the view was
    zoomed asymmetrically - with its label half off the plot. Centred, the
    point and the label naming it sit together in the middle.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to set the limits on
    moment : datetime.datetime
        The instant to centre on
    span_seconds : float
        Total width of view to show around it
    """

    half = timedelta(seconds=max(span_seconds, 1) / 2.0)
    ax.set_xlim(date2num(moment - half), date2num(moment + half))

    return None


# rotations tried, in order, when the boxplot's category labels don't fit
# horizontally. Each is anchored at the tick (ha="right") so the label reads
# with its right end at the tick and its left end lower, the conventional
# way a rotated x-tick label is drawn
_BOXPLOT_LABEL_ROTATIONS = (0, 15, 30, 45, 60, 75, 90)

# rotation used when the labels are switched on by hand despite none of the
# above fitting - a tilt rather than the near-vertical top of that list, so
# a long label is spread over some width instead of clipping straight down
_BOXPLOT_LABEL_FALLBACK_ROTATION = 45


def _boxplot_labels_fit(ax, renderer, container):
    """
    Whether every x-tick label currently installed on `ax` sits fully inside
    `container` and doesn't overlap any of its neighbours.

    Real, already-rendered bounding boxes are measured rather than predicted
    from font metrics - see _timeseries_labels_fit() for why a predicted
    width is not trusted here.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis whose x-tick labels are checked.
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure the labels against.
    container : matplotlib.transforms.Bbox
        Box every label must stay inside - only its horizontal extent
        generally matters, as a label is drawn below its own axis in any
        case, so the box is usually bounded to the axis's own column but
        left open vertically.

    Returns
    -------
    bool
        Whether they all fit.
    """

    boxes = [
        label.get_window_extent(renderer) for label in ax.xaxis.get_majorticklabels()
    ]

    for box in boxes:
        if (box.x0 < container.x0) or (box.x1 > container.x1):
            return False
        if (box.y0 < container.y0) or (box.y1 > container.y1):
            return False

    for i, box in enumerate(boxes):
        for other in boxes[i + 1 :]:
            if box.overlaps(other):
                return False

    return True


def fit_boxplot_xticklabels(
    ax, xticks, xtick_labels, xtick_params, xticklabel_params, forced=None
):
    """
    Show the boxplot's category labels along the x-axis if they can be made
    to fit, and hide them altogether otherwise.

    Horizontal is tried first, then progressively steeper rotations (see
    _BOXPLOT_LABEL_ROTATIONS), stopping at the first one under which no
    label runs off the figure and no two neighbouring labels overlap. If
    none of them fit, the labels are hidden rather than left to clash or
    spill off screen - the dashboard's panels are narrow enough, and vary
    enough with the layout, for that to happen often.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to set the boxplot's x-ticks/labels on.
    xticks : list
        Tick positions.
    xtick_labels : list
        Label text for each tick.
    xtick_params : dict
        Base tick parameters, as passed to xaxis.set_tick_params() -
        "rotation" and "labelbottom" are overwritten per candidate tried.
    xticklabel_params : dict
        Base label parameters, as passed to set_xticklabels() - "ha" and
        "rotation_mode" are overwritten per candidate tried.
    forced : bool, optional
        Skip the fitting search and show (True) or hide (False) the labels
        regardless of whether they fit - a manual override of what the
        automatic search would otherwise decide. Forcing them on when
        nothing actually fits falls back to the steepest rotation tried, as
        the most compact of the lot. Left as None, the fit is decided
        automatically.

    Returns
    -------
    bool
        Whether the labels ended up shown.
    """

    def install(rotation):
        params = copy.deepcopy(xtick_params)
        label_params = copy.deepcopy(xticklabel_params)
        params["rotation"] = rotation
        params["labelbottom"] = True
        if rotation == 0:
            label_params["ha"] = "center"
            label_params.pop("rotation_mode", None)
        else:
            label_params["ha"] = "right"
            label_params["rotation_mode"] = "anchor"
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels, rotation=rotation, **label_params)
        ax.xaxis.set_tick_params(**params)

    def hide():
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)
        ax.xaxis.set_tick_params(labelbottom=False)

    if forced is False:
        hide()
        return False

    # a single label cannot clash with anything but itself, and reads best
    # horizontal
    if len(xtick_labels) <= 1:
        install(0)
        return True

    # one real draw, to settle the figure's own layout - reused as the
    # renderer for every candidate rotation tried below (see
    # _timeseries_labels_fit() for why this measures real, rendered boxes
    # rather than predicting them)
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()

    # the dashboard packs several independent panels onto one shared figure,
    # so a label is "off screen" once it spills past this axis's own column
    # - into a neighbouring panel, or past the window's edge if it is the
    # outermost one - not merely once it clears the whole figure, which
    # would let it run straight through whatever panel sits next to it.
    # Horizontally bounded to the axis itself; vertically bounded to the
    # whole figure, as a label is drawn below its own axis in any case but
    # a steep rotation on a long label can still run past the bottom of a
    # dashboard panel tucked into the last row
    axis_box = ax.get_window_extent(renderer)
    container = mtransforms.Bbox.from_extents(
        axis_box.x0, ax.figure.bbox.y0, axis_box.x1, ax.figure.bbox.y1
    )

    for rotation in _BOXPLOT_LABEL_ROTATIONS:
        install(rotation)
        ax.figure.canvas.draw()
        if _boxplot_labels_fit(ax, renderer, container):
            return True

    # nothing fitted whole - forcing them on anyway falls back to a
    # moderate tilt rather than the steepest rotation just tried: that runs
    # close to vertical, which for a long label clips far more of it against
    # the bottom of the panel than a 45-degree tilt does, for no gain in
    # width. Not the most compact option, but the most legible of the ones
    # that clip
    if forced:
        install(_BOXPLOT_LABEL_FALLBACK_ROTATION)
        return True

    hide()
    return False


def compute_timeseries_xticks(
    ax, left, right, automatic_max_ticks=6, min_gap_pixels=12, data_start=None,
    data_end=None, data_resolution_seconds=None, n_ticks=None, force_edge_ticks=False,
):
    """
    Set x-axis tick positions and labels for a timeseries date range directly
    on `ax`, evenly spaced and landing on calendar boundaries that mean
    something - a month start, a midnight, an hour on the clock - rather than
    wherever a generic locator would put them.

    Every label shown is one step of a single stride from its neighbours, and
    a stride is taken whole or not at all: no tick is ever dropped to make
    the rest of a finer one fit. An earlier version thinned a tier that
    almost fitted, which left runs like twelve months with one missing and a
    double-width gap where it had been. Sparser but even reads far better
    than denser but gap-toothed, so the stride simply widens until the whole
    of it fits.

    Where the ticks land follows from that: the start and end of the range
    are labelled whenever they fall on the stride, which they do for the
    ranges usually loaded - a whole number of months, or days counted from
    the first - and are left unlabelled when the view has been zoomed
    somewhere that does not line up. Nothing is forced onto an edge, as a
    forced edge label sits at its own distance from its neighbour and is
    exactly the unevenness this is here to avoid. Day strides count from the
    first day of the loaded data so the data's own start and end fall on the
    grid; hours count from the clock, so a sub-daily view lands on 12h, 6h,
    3h and so on, always including midnight.

    Density adapts to the axis's actual, current pixel width rather than a
    fixed guess at how many ticks "should" fit: strides are tried from finest
    to coarsest, and each is really installed, drawn and measured to see
    whether all of its labels clear one another. So a wide panel naturally
    ends up with more, closer-together ticks than a narrow one showing the
    same span, without either ever clashing.

    Three earlier versions of the fitting step measured candidates *before*
    they were on screen - a fraction of the total data range, then a
    prediction of each label's rendered pixel width - and still let labels
    clash on a real machine despite passing the same kind of check in a test
    environment: a predicted width is only as good as the font/DPI
    assumptions behind it, and those can differ (font substitution, HiDPI
    scaling, hinting) between wherever this gets tested and where it runs. So
    this measures real, already-rendered label bounding boxes instead.

    The view is widened a little at each end where the outermost labels need
    it, so they sit centred under their ticks rather than being shunted
    sideways to stay on the axis - which read as one off-centre label at the
    end of an otherwise even row.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The real axis to set ticks on - candidates are installed on it and it
        is actually drawn so each label's real rendered bounding box can be
        read back.
    left : datetime.datetime
        Start of the visible x-axis range.
    right : datetime.datetime
        End of the visible x-axis range.
    automatic_max_ticks : int, default 6
        Ceiling on how many ticks a stride may offer to even be considered.
        Not a target: the stride actually used is whichever fits this ceiling
        *and* the axis's real available width, so the result is commonly
        fewer on a narrow panel and can be noticeably more on a wide one.
        Ignored when `n_ticks` is set.
    min_gap_pixels : float, default 12
        Minimum gap, in pixels, required between two adjacent labels' edges.
    data_start, data_end : datetime.datetime, optional
        The true, configured start/end of the loaded data (not the current
        view). data_start is what multi-day strides are counted from, so the
        data's own start and end land on the grid.
    data_resolution_seconds : float, optional
        The loaded data's own sampling interval, in seconds (e.g. 3600 for
        hourly). When the current view is zoomed to no wider than this, there
        is at most one real observation actually in it, so a single tick at
        that observation's own time is shown instead. Left as None, this
        collapse never happens.
    n_ticks : int, optional
        Manual tick count - evenly spaced on whole hours between the visible
        data's first and last hour, instead of the calendar-aligned automatic
        search, and `automatic_max_ticks` is ignored. None (the default) keeps
        the automatic behaviour.
    force_edge_ticks : bool, default False
        Always label the visible data's first and last whole hour, even where
        they fall off whatever stride was chosen (the automatic search otherwise
        leaves an edge unlabelled rather than space it unevenly - see above).

    Returns
    -------
    xticks : list of datetime.datetime
        The tick positions actually set on `ax`, in order.
    """

    full_precision = "%Y-%m-%d %H:%M:%S"

    if left >= right:
        dates = [left, right] if left < right else [left]
        ax.xaxis.set_ticks(dates, labels=[d.strftime(full_precision) for d in dates])
        if len(dates) == 1:
            # nothing but a single instant to show, so it goes in the middle
            # with a sample interval either side of it - a view of no width
            # is left to matplotlib to expand however it sees fit
            _centre_view_on(ax, dates[0], (data_resolution_seconds or 3600) * 2)
        return dates

    # the plotted data's own start/end where the view still shows them, not
    # the margin-padded view limits - rounded inwards to whole hours, as no
    # finer resolution is ever used
    plotted_extent = _plotted_time_extent(ax)
    if plotted_extent is not None:
        edge_left = max(left, plotted_extent[0])
        edge_right = min(right, plotted_extent[1])
    else:
        edge_left = data_start if (data_start and left <= data_start <= right) else left
        edge_right = data_end if (data_end and left <= data_end <= right) else right
    edge_left = pd.Timestamp(edge_left).ceil("h").to_pydatetime()
    edge_right = pd.Timestamp(edge_right).floor("h").to_pydatetime()
    have_hour_edges = edge_left <= edge_right

    if (n_ticks is not None) and have_hour_edges:
        dates = sorted(
            {
                d.round("h").to_pydatetime()
                for d in pd.date_range(edge_left, edge_right, periods=max(n_ticks, 2))
            }
        )
        xlim = ax.get_xlim()
        ax.xaxis.set_ticks(dates, labels=_format_free_timeseries_ticks(dates))
        ax.set_xlim(*xlim)
        _set_timeseries_tick_alignment(ax.xaxis.get_majorticklabels())
        return dates

    if (
        data_resolution_seconds
        and (right - left).total_seconds() <= data_resolution_seconds * 1.01
    ):
        # the view is no wider than a single sample interval, so show that
        # one observation's own instant rather than two edge labels either
        # side of a gap with nothing in it. The tick goes on the sample
        # boundary nearest the middle of the view, not the view's midpoint,
        # as the label has to sit under the plotted point it names. Only
        # boundaries inside the view are eligible, as matplotlib silently
        # widens the axis to fit a tick outside it
        centre = left + (right - left) / 2
        epoch = datetime(left.year, 1, 1)
        step = data_resolution_seconds
        steps_before = int((centre - epoch).total_seconds() // step)
        in_view = [
            boundary
            for boundary in (
                epoch + timedelta(seconds=steps_before * step),
                epoch + timedelta(seconds=(steps_before + 1) * step),
            )
            if left <= boundary <= right
        ]
        if in_view:
            snapped = min(
                in_view, key=lambda boundary: abs((boundary - centre).total_seconds())
            )
            ax.xaxis.set_ticks(
                [snapped], labels=[_format_single_observation_tick(snapped)]
            )
            # the one observation goes in the middle of the panel, so it and
            # the label naming it are read together rather than sitting off
            # to whichever side the view happened to be zoomed to
            _centre_view_on(ax, snapped, (right - left).total_seconds())
            return [snapped]

    # one real draw, to settle the axis's own layout - reused as the
    # renderer for every candidate set tried below instead of redrawing
    # the whole figure (data, other axes, everything) each time; see
    # _timeseries_labels_fit() for why that's still an accurate measure
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()

    # the view's numeric x-limits, re-applied after every set_ticks() call
    # below - matplotlib silently expands the axis's data limits to fit any
    # tick outside the current view, even with autoscale off
    xlim = ax.get_xlim()

    # prime the axis's tick-label object pool before the search takes its
    # first measurement: the first Text objects matplotlib hands out measure
    # a few pixels off from what every later call confirms, occasionally
    # enough to make the first tier tried look like it fits when it doesn't
    ax.xaxis.set_ticks([left, right], labels=["", ""])
    ax.set_xlim(*xlim)
    for label in ax.xaxis.get_majorticklabels():
        label.get_window_extent(renderer)

    def choose_ticks(current_xlim):
        """
        Pick the ticks to show, given the view's current limits.

        A tier is taken or left whole - no tick is ever dropped from one to
        make the rest fit, which is what left a run of months with one
        missing and the spacing visibly uneven. Sparser but even beats denser
        but gap-toothed.
        """

        grids = {}
        best_calendar = None
        best_calendar_count = 0
        best_plain = None
        best_plain_count = 0
        for step_kind, multiple in _TIMESERIES_TICK_STEPS:
            step_ticks = [
                dt
                for dt in _aligned_timeseries_ticks(
                    left, right, step_kind, multiple, anchor=data_start
                )
                if left <= dt <= right
            ]
            if len(step_ticks) < 2 or len(step_ticks) > automatic_max_ticks:
                continue

            texts = [_format_timeseries_tick(dt, step_kind) for dt in step_ticks]
            if not _timeseries_labels_fit(
                ax, renderer, step_ticks, texts, min_gap_pixels, current_xlim
            ):
                continue

            # a month start, or the 1st and 15th, are landmarks a reader
            # recognises, where a plain day stride lands on dates of no
            # significance - so the two are tracked apart and weighed up below
            if step_kind in _MONTH_DAY_GRIDS:
                grids[step_kind] = (step_ticks, texts, step_kind, multiple)
            elif step_kind in _CALENDAR_TICK_KINDS:
                if len(step_ticks) > best_calendar_count:
                    best_calendar = (step_ticks, texts, step_kind, multiple)
                    best_calendar_count = len(step_ticks)
            elif len(step_ticks) > best_plain_count:
                best_plain = (step_ticks, texts, step_kind, multiple)
                best_plain_count = len(step_ticks)

        # the coarsest grid that can fill the view wins, as the coarser it
        # is the more even it is: the 1st and the 15th first, then the 1st,
        # 10th and 20th once a view is too short to hold three of those,
        # then every fifth day. Below all of them - a week or less, where
        # even the finest grid puts only a label or two on the axis - a
        # plain stride counted from the data takes over, and above them one
        # label per month or per quarter does
        for grid_kind in _MONTH_DAY_GRIDS:
            grid = grids.get(grid_kind)
            if grid is not None and len(grid[0]) >= _MIN_MONTH_GRID_TICKS:
                return grid

        if best_calendar is not None:
            return best_calendar

        return best_plain

    # choosing the ticks and making room for the end labels pull against each
    # other - widening the view to fit a label spreads the same pixels over a
    # longer span, which pulls every tick closer to its neighbours and can
    # crowd labels that fitted a moment ago. So the two are run in turn until
    # the choice survives the room made for it
    kept_dates, kept_texts = None, None
    for _ in range(_EDGE_LABEL_FIT_PASSES):
        xlim = ax.get_xlim()
        chosen = choose_ticks(xlim)
        if chosen is None:
            break

        kept_dates, kept_texts = chosen[0], chosen[1]
        ax.xaxis.set_ticks(kept_dates, labels=kept_texts)
        ax.set_xlim(*xlim)
        _set_timeseries_tick_alignment(ax.xaxis.get_majorticklabels())
        _widen_for_edge_labels(ax, renderer, kept_dates)

        widened_xlim = ax.get_xlim()
        if widened_xlim == xlim:
            break
        if _timeseries_labels_fit(
            ax, renderer, kept_dates, kept_texts, min_gap_pixels, widened_xlim
        ):
            break

    if kept_dates is None:
        # nothing fitted whole, at any resolution - fall back to the two ends
        # of the view, which is the least that still says what is being shown
        kept_dates = [left, right]
        kept_texts = [left.strftime(full_precision), right.strftime(full_precision)]

    if force_edge_ticks and have_hour_edges:
        kept = [
            (date, text)
            for date, text in zip(kept_dates, kept_texts)
            if edge_left <= date <= edge_right
        ]
        kept_dates = [date for date, _ in kept]
        kept_texts = [text for _, text in kept]
        if not kept_dates or kept_dates[0] != edge_left:
            kept_dates.insert(0, edge_left)
            kept_texts.insert(0, _format_free_timeseries_ticks([edge_left])[0])
        if kept_dates[-1] != edge_right:
            kept_dates.append(edge_right)
            kept_texts.append(_format_free_timeseries_ticks([edge_right])[0])

    ax.xaxis.set_ticks(kept_dates, labels=kept_texts)
    ax.set_xlim(*ax.get_xlim())
    _set_timeseries_tick_alignment(ax.xaxis.get_majorticklabels())

    return kept_dates


def harmonise_xy_lims_paradigm(
    read_instance,
    canvas_instance,
    relevant_axs,
    base_plot_type,
    plot_characteristics,
    plot_options,
    xlim=None,
    ylim=None,
    relim=False,
    autoscale=False,
    autoscale_x=False,
    autoscale_y=False,
    bias_centre=False,
    harmonise=True,
):
    """
    Harmonise x and y axis limits across a set of axes for a given plot type,
    unless axis limits have been manually defined.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    canvas_instance : object
        Instance of class Canvas or Report.
    relevant_axs : list, dict or object
        Axes to harmonise limits for. Can be a single axis, a list of axes,
        or a dict mapping temporal resolutions to axes.
    base_plot_type : str
        Plot type, without statistical overlays (e.g., 'scatter', 'timeseries', 'periodic').
    plot_characteristics : dict
        Plot characteristics, including optional 'xlim', 'ylim', 'equal_aspect',
        'xtick_alteration', and 'margin_padding'.
    plot_options : list
        Active plot options (e.g., ['bias', ...]).
    xlim : dict or None
        Optional x-axis limits to set.
    ylim : dict or None
        Optional y-axis limits to set.
    relim : bool, default False
        If True, recompute axis limits from data before harmonisation.
    autoscale : bool, default False
        If True, autoscale both x and y axes to the data.
    autoscale_x : bool, default False
        If True, autoscale x-axis to the data.
    autoscale_y : bool, default False
        If True, autoscale y-axis to the data.
    bias_centre : bool, default False
        If True and 'bias' in plot_options, centre y-axis limits at zero.
    harmonise : bool, default True
        If True, harmonise axes across the paradigm.
    """

    # periodic_relevant_temporal_resolutions is set once data has been
    # resampled (see resample() in statistics.py), which can lag behind the
    # first time a periodic axis is interacted with - guard rather than crash
    if base_plot_type in ["periodic", "periodic-violin"] and not hasattr(
        read_instance, "periodic_relevant_temporal_resolutions"
    ):
        active_resolution = getattr(read_instance, "active_resolution", None) or getattr(
            read_instance, "resolution", None
        )
        if active_resolution:
            read_instance.periodic_relevant_temporal_resolutions = (
                get_periodic_relevant_temporal_resolutions(active_resolution)
            )

    # initialise arrays to save lower and upper limits in all axes
    all_xlim_lower = []
    all_xlim_upper = []
    all_ylim_lower = []
    all_ylim_upper = []

    # initialise variables for setting axis limits
    xlim_min = None
    xlim_max = None
    ylim_min = None
    ylim_max = None

    # transform axis dict or str to list. A dict's own keys drive the
    # periodic resolution mapping below, so a caller can pass just the one
    # sub-axis a zoom changed rather than have harmonisation pull in the
    # other panels' limits (see harmonise_changed_axis() in toolbar.py)
    dict_resolutions = None
    if not isinstance(relevant_axs, list):
        # if changes only apply to one axis, put it in list
        if not isinstance(relevant_axs, dict):
            relevant_axs = [relevant_axs]
        # transform dictionaries into lists
        else:
            dict_resolutions = list(relevant_axs.keys())
            relevant_axs = [relevant_axs[k] for k in dict_resolutions]

    # get mapped resolution per axis for periodic plots
    if base_plot_type in ["periodic", "periodic-violin"]:
        if dict_resolutions is not None:
            mapped_resolutions = dict_resolutions
        else:
            mapped_resolutions = read_instance.periodic_relevant_temporal_resolutions * (
                int(
                    len(relevant_axs)
                    / len(read_instance.periodic_relevant_temporal_resolutions)
                )
            )

    # remove any axes from relevant_axs which are not active (only for report and library),
    # and any that are None - a caller can end up passing one (e.g. a
    # stale/mismatched axis-to-plot-type lookup at the toolbar layer)
    # and there's nothing to harmonise for it anyway
    if read_instance.mode in ["report", "library"]:
        relevant_axs_active = []
        mapped_resolutions_active = []
        for ax_ii, ax in enumerate(relevant_axs):
            if ax is not None and ax.axison:
                relevant_axs_active.append(ax)
                if base_plot_type in ["periodic", "periodic-violin"]:
                    mapped_resolutions_active.append(mapped_resolutions[ax_ii])
    else:
        if base_plot_type in ["periodic", "periodic-violin"]:
            relevant_axs_active = []
            mapped_resolutions_active = []
            for ax_ii, ax in enumerate(relevant_axs):
                if ax is not None:
                    relevant_axs_active.append(ax)
                    mapped_resolutions_active.append(mapped_resolutions[ax_ii])
        else:
            relevant_axs_active = [ax for ax in relevant_axs if ax is not None]

    # get lower and upper limits across all relevant axes
    for ax in relevant_axs_active:
        if "equal_aspect" in plot_characteristics:
            if plot_characteristics["equal_aspect"]:
                set_equal_axes(ax, plot_options, plot_characteristics, base_plot_type)
        else:
            ax.set_aspect("auto")

        if relim:
            ax.relim(visible_only=True)
        if autoscale:
            ax.autoscale(tight=False)
        if autoscale_x:
            ax.autoscale(axis="x", tight=False)
        if autoscale_y:
            ax.autoscale(axis="y", tight=False)

        if (xlim is None) and ("xlim" not in plot_characteristics):
            if base_plot_type not in [
                "timeseries",
                "boxplot",
                "periodic",
                "periodic-violin",
            ]:
                xlim_lower, xlim_upper = ax.get_xlim()
            elif base_plot_type in ["timeseries", "scatter"]:
                if base_plot_type == "timeseries":
                    # the axis's actual current view, not a margin-stripped
                    # derivative - get_no_margin_lim() assumes the configured
                    # margin is baked into xlim, which is only true right
                    # after an autoscale. An explicit zoom/pan sets xlim with
                    # no margin, so subtracting one back out shrinks the
                    # visible range and mislabels the start/end ticks
                    xlim_lower, xlim_upper = ax.get_xlim()
                else:
                    xlim_lower, xlim_upper = get_no_margin_lim(ax, "xlim")
                try:
                    xlim_lower = num2date(xlim_lower).replace(tzinfo=None)
                    xlim_upper = num2date(xlim_upper).replace(tzinfo=None)
                except ValueError:
                    continue

            if base_plot_type not in ["boxplot", "periodic", "periodic-violin"]:
                all_xlim_lower.append(xlim_lower)
                all_xlim_upper.append(xlim_upper)

        if (ylim is None) and ("ylim" not in plot_characteristics):
            ylim_lower, ylim_upper = ax.get_ylim()
            all_ylim_lower.append(ylim_lower)
            all_ylim_upper.append(ylim_upper)

    # get minimum and maximum from all axes and set limits
    for ax_ii, ax in enumerate(relevant_axs_active):
        # get xlim
        if ax_ii == 0:
            set_xlim = False
            if (
                (xlim is None)
                and ("xlim" not in plot_characteristics)
                and (len(all_xlim_lower) > 0)
                and (len(all_xlim_upper) > 0)
            ):
                if base_plot_type not in ["boxplot", "periodic", "periodic-violin"]:
                    xlim_min = np.min(all_xlim_lower)
                    xlim_max = np.max(all_xlim_upper)
                    xlim = {"left": xlim_min, "right": xlim_max}
                    if harmonise:
                        set_xlim = True
            elif "xlim" in plot_characteristics:
                xlim = plot_characteristics["xlim"]
                set_xlim = True

        # set xlim
        if set_xlim and (
            base_plot_type
            not in ["timeseries", "boxplot", "periodic", "periodic-violin"]
        ):
            if isinstance(xlim, dict):
                ax.set_xlim(**xlim)
            else:
                ax.set_xlim(xlim)

        # get ylim
        if ax_ii == 0:
            set_ylim = False
            if (
                (ylim is None)
                and ("ylim" not in plot_characteristics)
                and (len(all_ylim_lower) > 0)
                and (len(all_ylim_upper) > 0)
            ):
                ylim_min = np.min(all_ylim_lower)
                ylim_max = np.max(all_ylim_upper)
                # if have bias_centre option, centre around zero
                if ("bias" in plot_options) & (bias_centre):
                    if np.abs(np.max(all_ylim_upper)) >= np.abs(np.min(all_ylim_lower)):
                        ylim_min = -np.abs(np.max(all_ylim_upper))
                        ylim_max = np.abs(np.max(all_ylim_upper))
                    elif np.abs(np.max(all_ylim_upper)) < np.abs(
                        np.min(all_ylim_lower)
                    ):
                        ylim_min = -np.abs(np.min(all_ylim_lower))
                        ylim_max = np.abs(np.min(all_ylim_lower))
                ylim = {"bottom": ylim_min, "top": ylim_max}
                if harmonise:
                    set_ylim = True
            elif "ylim" in plot_characteristics:
                ylim = plot_characteristics["ylim"]
                set_ylim = True

        # set ylim
        if set_ylim:
            if isinstance(ylim, dict):
                ax.set_ylim(**ylim)
            else:
                ax.set_ylim(ylim)

    # get minimum and maximum from all axes and set limits for periodic plots
    if base_plot_type in ["periodic", "periodic-violin"]:
        if (xlim is None) and ("xlim" not in plot_characteristics):
            for temporal_resolution, sub_ax in zip(
                mapped_resolutions_active, relevant_axs_active
            ):
                # adjust plot x axis to have correct margin on edges
                xlim_lower, xlim_upper = sub_ax.get_xlim()
                first_valid_x = canvas_instance.periodic_xticks[temporal_resolution][
                    (
                        np.abs(
                            canvas_instance.periodic_xticks[temporal_resolution]
                            - xlim_lower
                        )
                    ).argmin()
                ]
                last_valid_x = canvas_instance.periodic_xticks[temporal_resolution][
                    (
                        np.abs(
                            canvas_instance.periodic_xticks[temporal_resolution]
                            - xlim_upper
                        )
                    ).argmin()
                ]
                if temporal_resolution == "hour":
                    xlim_lower = first_valid_x - 0.65
                    xlim_upper = last_valid_x + 0.65

                    # the "hour" axis is normally shown at every 3rd hour
                    # (see format_axis()), a fixed step set once when the
                    # plot is built - zooming into a handful of hours could
                    # land between two of those and show no ticks at all, so
                    # recompute the step from the hours actually in view
                    visible_hours = last_valid_x - first_valid_x + 1
                    if visible_hours <= 8:
                        hour_step = 1
                    elif visible_hours <= 16:
                        hour_step = 2
                    else:
                        hour_step = 3
                    sub_ax.set_xticks(
                        canvas_instance.periodic_xticks[temporal_resolution][
                            ::hour_step
                        ]
                    )
                elif temporal_resolution == "dayofweek":
                    xlim_lower = first_valid_x - 0.55
                    xlim_upper = last_valid_x + 0.55
                elif temporal_resolution == "month":
                    xlim_lower = first_valid_x - 0.55
                    xlim_upper = last_valid_x + 0.55
                xlim = {"left": xlim_lower, "right": xlim_upper}
                sub_ax.set_xlim(**xlim)
        elif "xlim" in plot_characteristics:
            xlim = plot_characteristics["xlim"]
            for temporal_resolution, sub_ax in zip(
                mapped_resolutions_active, relevant_axs_active
            ):
                if isinstance(xlim, dict):
                    sub_ax.set_xlim(**xlim)
                else:
                    sub_ax.set_xlim(xlim)

        # if harmonisation is off, and ylim not manually set,
        # ensure harmonisation is at least done for a plot across resolutions
        if (not harmonise) and (not set_ylim):
            current_resolutions = []
            current_axs = []
            current_ylim_lower = []
            current_ylim_upper = []

            for ax_ii, (temporal_resolution, sub_ax) in enumerate(
                zip(mapped_resolutions_active, relevant_axs_active)
            ):
                # get temporal resolution of next axis
                if ax_ii != (len(relevant_axs_active) - 1):
                    next_temporal_resolution = mapped_resolutions_active[ax_ii + 1]
                else:
                    next_temporal_resolution = None

                # if resolution not yet in current resolutions, then add information for it
                if temporal_resolution not in current_resolutions:
                    ylim_lower, ylim_upper = sub_ax.get_ylim()
                    current_resolutions.append(temporal_resolution)
                    current_axs.append(sub_ax)
                    current_ylim_lower.append(ylim_lower)
                    current_ylim_upper.append(ylim_upper)

                # if next resolution already in current resolutions or on last axis, then set ylim for relevant axes
                if (next_temporal_resolution in current_resolutions) or (
                    ax_ii == (len(relevant_axs_active) - 1)
                ):
                    ylim_min = np.min(current_ylim_lower)
                    ylim_max = np.max(current_ylim_upper)
                    # if have bias_centre option, centre around zero
                    if ("bias" in plot_options) & (bias_centre):
                        if np.abs(np.max(current_ylim_upper)) >= np.abs(
                            np.min(current_ylim_lower)
                        ):
                            ylim_min = -np.abs(np.max(current_ylim_upper))
                            ylim_max = np.abs(np.max(current_ylim_upper))
                        elif np.abs(np.max(current_ylim_upper)) < np.abs(
                            np.min(current_ylim_lower)
                        ):
                            ylim_min = -np.abs(np.min(current_ylim_lower))
                            ylim_max = np.abs(np.min(current_ylim_lower))
                    ylim = {"bottom": ylim_min, "top": ylim_max}
                    for current_ax in current_axs:
                        current_ax.set_ylim(**ylim)

                    # reset lists
                    current_resolutions = []
                    current_axs = []
                    current_ylim_lower = []
                    current_ylim_upper = []

    # get minimum and maximum from all axes and set limits for timeseries
    elif base_plot_type == "timeseries":
        if (plot_characteristics["xtick_alteration"]["define"]) and (xlim):
            # get left and right
            if isinstance(xlim, dict):
                left = xlim["left"]
                right = xlim["right"]
            else:
                left = xlim[0]
                right = xlim[1]

            if left == right:
                # a degenerate (zero-width) range has no meaningful ticks to
                # compute - let matplotlib pick, as in the unresolved xlim
                # case below
                for ax in relevant_axs_active:
                    ax.xaxis.set_major_locator(mpl.dates.AutoDateLocator())
                    ax.xaxis.set_major_formatter(
                        mpl.dates.ConciseDateFormatter(ax.xaxis.get_major_locator())
                    )
                return

            if read_instance.daily_forecast:
                # forecast-day-numbered labels ("Day1 3h") need day-aligned
                # boundaries to count from, so round to the nearest day
                if left.hour >= 12:
                    forecast_left = datetime(
                        left.year, left.month, left.day
                    ) + timedelta(days=1)
                else:
                    forecast_left = datetime(left.year, left.month, left.day)
                if right.hour >= 12:
                    forecast_right = datetime(
                        right.year, right.month, right.day
                    ) + timedelta(days=1)
                else:
                    forecast_right = datetime(right.year, right.month, right.day)

                forecast_hours = (
                    forecast_right - forecast_left
                ).total_seconds() / 3600
                if forecast_hours <= 24:
                    freq = "3h"
                elif forecast_hours <= 48:
                    freq = "6h"
                else:
                    freq = "12h"
                xticks = pd.date_range(forecast_left, forecast_right, freq=freq)
                xticklabels = []
                start_pd_dt = xticks[0]
                for pd_dt in xticks:
                    pd_dt_diff = pd_dt - start_pd_dt
                    day = (
                        pd_dt_diff.days
                        + 1
                        + (read_instance.active_forecast_days[0] - 1)
                    )
                    hour = pd_dt.strftime("%H")
                    if int(hour) == 0:
                        label = "Day{} {}h".format(day, hour)
                    else:
                        label = "{}h".format(hour)
                    xticklabels.append(label)
                for ax in relevant_axs_active:
                    ax.xaxis.set_ticks(xticks, labels=xticklabels)

            else:
                # ticks snap to calendar-aligned resolutions (see
                # compute_timeseries_xticks()) unless "n_ticks" overrides that
                xtick_alteration = plot_characteristics["xtick_alteration"]
                automatic_max_ticks = xtick_alteration["automatic_max_ticks"]
                n_ticks = xtick_alteration.get("n_ticks")
                force_edge_ticks = xtick_alteration.get("force_edge_ticks", False)
                data_start = _parse_yyyymmdd(getattr(read_instance, "start_date", None))
                data_end = _parse_yyyymmdd(getattr(read_instance, "end_date", None))
                active_resolution = getattr(read_instance, "active_resolution", None) or getattr(
                    read_instance, "resolution", None
                )
                data_resolution_seconds = _RESOLUTION_SECONDS.get(active_resolution)
                for ax in relevant_axs_active:
                    compute_timeseries_xticks(
                        ax, left, right, automatic_max_ticks=automatic_max_ticks,
                        data_start=data_start, data_end=data_end,
                        data_resolution_seconds=data_resolution_seconds,
                        n_ticks=n_ticks, force_edge_ticks=force_edge_ticks,
                    )

            # pad the margins
            for ax in relevant_axs_active:
                ax.margins(**plot_characteristics["margin_padding"])

        elif plot_characteristics["xtick_alteration"]["define"]:
            # xlim couldn't be resolved for any axis this call - don't leave
            # ticks as whatever a previous, differently zoomed call set:
            # those positions are unlikely to fall inside a narrowed view,
            # showing no ticks at all rather than just imprecise ones
            for ax in relevant_axs_active:
                ax.xaxis.set_major_locator(mpl.dates.AutoDateLocator())
                ax.xaxis.set_major_formatter(
                    mpl.dates.ConciseDateFormatter(ax.xaxis.get_major_locator())
                )


def set_axis_title(read_instance, relevant_axis, title, plot_characteristics):
    """
    Set the title of a plot axis.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    relevant_axis : object, list or dict
        Axis (or axes) to set the title for. Can be a single axis, a list of axes,
        or a dict mapping temporal resolutions to axes.
    title : str
        Title to set on the axis.
    plot_characteristics : dict
        Plot characteristics containing formatting options for the axis title.
    """

    # return if title is empty str
    if title == "":
        return

    # get appropriate axis for plotting label for plots with multiple sub-axes
    axs_to_set_title = []
    if isinstance(relevant_axis, dict):
        # reorder dict to show axis title in monthly plot and not in DoW for daily plots
        relevant_dict = {
            key: relevant_axis[key] for key in ["hour", "month", "dayofweek"]
        }
        for relevant_temporal_resolution, sub_ax in relevant_dict.items():
            if (
                relevant_temporal_resolution
                in read_instance.periodic_relevant_temporal_resolutions
            ):
                axs_to_set_title.append(sub_ax)
                break
    elif isinstance(relevant_axis, list):
        axs_to_set_title.append(relevant_axis[0])
    else:
        axs_to_set_title.append(relevant_axis)

    # set title for appropriate axes
    axis_title_characteristics = copy.deepcopy(plot_characteristics["axis_title"])
    axis_title_characteristics["label"] = title
    for relevant_axis in axs_to_set_title:
        relevant_axis.set_title(**axis_title_characteristics)


def set_axis_label(
    relevant_axis,
    label_ax,
    label,
    plot_characteristics,
    relevant_temporal_resolutions=None,
):
    """
    Set the label of a plot axis.

    Parameters
    ----------
    relevant_axis : object, list or dict
        Axis (or axes) to set the label for. Can be a single axis or a dict mapping
        temporal resolutions to axes.
    label_ax : str
        Which axis to set the label for: 'x' or 'y'.
    label : str
        Label text to set.
    plot_characteristics : dict
        Plot characteristics containing formatting options under 'xlabel' or 'ylabel'.
    relevant_temporal_resolutions : list, optional
        List of temporal resolutions to include when setting labels on dict axes.
    """

    # return if label is empty str
    if label == "":
        return

    # define default argument mutables
    if relevant_temporal_resolutions is None:
        relevant_temporal_resolutions = ["hour", "month"]

    # get appropriate axis for plotting label for plots with multiple sub-axes (hour and month axes)
    axs_to_set_label = []
    if isinstance(relevant_axis, dict):
        for relevant_temporal_resolution, sub_ax in relevant_axis.items():
            if relevant_temporal_resolution in relevant_temporal_resolutions:
                axs_to_set_label.append(sub_ax)
            # remove day of week axis label if setting ylabel
            if (relevant_temporal_resolution == "dayofweek") & (label_ax == "y"):
                sub_ax.yaxis.set_tick_params(which="both", labelleft=False)
                sub_ax.set_ylabel("")
    else:
        axs_to_set_label.append(relevant_axis)

    # set label for appropriate axes
    for relevant_axis in axs_to_set_label:
        if label_ax == "x":
            axis_label_characteristics = copy.deepcopy(plot_characteristics["xlabel"])
            axis_label_characteristics["xlabel"] = label
            relevant_axis.set_xlabel(**axis_label_characteristics)
        elif label_ax == "y":
            axis_label_characteristics = copy.deepcopy(plot_characteristics["ylabel"])
            axis_label_characteristics["ylabel"] = label
            relevant_axis.set_ylabel(**axis_label_characteristics)


def format_plot_options(
    read_instance,
    canvas_instance,
    relevant_axs,
    relevant_data_labels,
    networkspeci,
    base_plot_type,
    plot_type,
    plot_options,
    map_extent=False,
    chunk_stat=None,
    chunk_resolution=None,
):
    """
    Function that handles formatting of a plot axis,
    based on given plot options.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report containing data and configuration.
    canvas_instance : object
        Instance of class Canvas or Report providing plot characteristics and layout.
    relevant_axs : list, dict or object
        Axes to apply formatting to. Can be a single axis, a list of axes,
        or a dictionary of sub-axes for periodic plots.
    relevant_data_labels : list
        Data labels corresponding to plotted data for each axis.
    networkspeci : str
        Current networkspeci identifier.
    base_plot_type : str
        Base plot type without statistical information.
    plot_type : str
        Specific plot type with optional statistical information.
    plot_options : list
        List of plot options to apply, such as 'logx', 'logy', 'domain', 'annotate',
        'regression', 'smooth', 'threshold'.
    map_extent : list, optional
        Spatial bounds for map plots [lon_min, lon_max, lat_min, lat_max].
    chunk_stat : str, optional
        Chunk statistic to use for smoothing or other statistical plot enhancements.
    chunk_resolution : str, optional
        Chunk resolution (e.g., 'day', 'month') for statistical calculations.
    """

    # transform axis dict or str to list
    if not isinstance(relevant_axs, list):
        # if changes only apply to one axis, put it in list
        if not isinstance(relevant_axs, dict):
            relevant_axs = [relevant_axs]
        # transform dictionaries into lists
        else:
            relevant_axs = [
                relevant_axs[relevant_temporal_resolution]
                for relevant_temporal_resolution in read_instance.periodic_relevant_temporal_resolutions
            ]
            relevant_data_labels = copy.deepcopy(relevant_data_labels) * len(
                read_instance.periodic_relevant_temporal_resolutions
            )

    # get zstat info (if any)
    (
        zstat,
        base_zstat,
        z_statistic_type,
        z_statistic_sign,
        z_statistic_period,
    ) = get_z_statistic_info(plot_type=plot_type)

    for relevant_ax_ii, relevant_ax in enumerate(relevant_axs):
        # log axes?
        if "logx" in plot_options:
            log_valid = log_validity(relevant_ax, "logx")
            if log_valid:
                log_axes(
                    relevant_ax, "logx", canvas_instance.plot_characteristics[plot_type]
                )
            else:
                msg = "It is not possible to log the x-axis "
                msg += "in {0} with negative values.".format(plot_type)
                show_message(read_instance, msg)

        if "logy" in plot_options:
            log_valid = log_validity(relevant_ax, "logy")
            if log_valid:
                log_axes(
                    relevant_ax, "logy", canvas_instance.plot_characteristics[plot_type]
                )
            else:
                msg = "It is not possible to log the y-axis "
                msg += "in {0} with negative values.".format(plot_type)
                show_message(read_instance, msg)

        # domain
        if "domain" in plot_options:
            if len(read_instance.data_labels) == 1:
                if (
                    read_instance.data_labels[0]
                    == read_instance.observations_data_label
                ):
                    msg = "'domain' plot option cannot be made as have no models."
                    show_message(read_instance, msg)
                    return
            model_domain(
                canvas_instance,
                relevant_ax,
                relevant_data_labels[relevant_ax_ii],
                map_extent,
            )

        # annotation
        if "annotate" in plot_options:
            if base_plot_type not in ["heatmap"]:
                annotation(
                    read_instance,
                    canvas_instance,
                    relevant_ax,
                    networkspeci,
                    relevant_data_labels[relevant_ax_ii],
                    base_plot_type,
                    canvas_instance.plot_characteristics[plot_type],
                    plot_options,
                    plot_z_statistic_sign=z_statistic_sign,
                )
                # annotate on first axis
                if base_plot_type in ["periodic", "periodic-violin"]:
                    break

        # regression line
        if "regression" in plot_options:
            linear_regression(
                read_instance,
                canvas_instance,
                relevant_ax,
                networkspeci,
                relevant_data_labels[relevant_ax_ii],
                base_plot_type,
                canvas_instance.plot_characteristics[plot_type],
                plot_options,
            )

        # smooth line
        if "smooth" in plot_options:
            smooth(
                read_instance,
                canvas_instance,
                relevant_ax,
                networkspeci,
                relevant_data_labels[relevant_ax_ii],
                base_plot_type,
                canvas_instance.plot_characteristics[plot_type],
                plot_options,
                chunk_stat,
                chunk_resolution,
            )

        # threshold line
        if "threshold" in plot_options:
            threshold(
                read_instance,
                canvas_instance,
                relevant_ax,
                networkspeci,
                base_plot_type,
                canvas_instance.plot_characteristics[plot_type],
            )


def map_feature_ink(map_template, feature="borders"):
    """
    A line colour for the map's borders and gridlines that stays visible
    against whatever land and ocean colours are currently set.

    Both are drawn over the basemap, so a single fixed colour cannot work
    for both a light and a dark one: the configured dark grey vanished
    entirely once the land and ocean were switched to their dark options.
    This picks a dark ink over a light basemap and a light ink over a dark
    one, judged on the perceived (luminance-weighted) brightness of the
    two basemap colours together rather than either alone, since borders
    and gridlines cross both.

    Parameters
    ----------
    map_template : dict
        The map's plot_characteristics_templates entry, holding the
        current land_polygon/ocean_polygon face colours.
    feature : {"borders", "gridlines"}, optional
        Which of the two to pick a colour for. Gridlines are a background
        reference the eye should be able to ignore, so they get a much
        softer tone than borders, which describe the map itself - drawing
        both at the same strength left the gridlines dominating the
        figure.

    Returns
    -------
    str
        Hex colour to draw the feature in.
    """

    def luminance(colour):
        # Rec. 709 relative luminance - matches how bright a colour
        # actually looks, unlike a plain mean of the channels
        red, green, blue = mpl.colors.to_rgb(colour)
        return 0.2126 * red + 0.7152 * green + 0.0722 * blue

    land, ocean = get_map_colours(map_template)
    try:
        brightness = (luminance(land) + luminance(ocean)) / 2
    except ValueError:
        # an unparseable colour shouldn't take the map down with it
        return "#8A8A8A" if feature == "gridlines" else "#4D4D4D"

    light_basemap = brightness > 0.5
    if feature == "gridlines":
        # only just enough separation from the basemap to be followed -
        # gridlines sit under everything else and are read by glancing at
        # them, never studied
        return "#9E9E9E" if light_basemap else "#6E7681"
    # not pure black/white at either end: full contrast makes borders
    # compete with the station data drawn on top of them, which is the
    # thing actually meant to stand out
    return "#3A3A3A" if light_basemap else "#C8C8C8"


def draw_map_features(canvas_instance, ax):
    """
    Add the map's ocean, land, and country border cartopy features to an
    axis, styled per canvas_instance.plot_characteristics_templates["map"].
    Only meaningful for the "providentia" map background (the default) - a
    custom background image or cartopy's shaded relief doesn't use these.

    Split out from format_axis() so it can be re-run on its own whenever the
    user changes land/ocean colour, border visibility, or map
    resolution from the map settings menu, without re-doing the rest of
    format_axis()'s one-time axis setup (which would duplicate gridlines).

    Parameters
    ----------
    canvas_instance : object
        Instance of class Canvas.
    ax : cartopy.mpl.geoaxes.GeoAxes
        Map axis to draw the features onto.

    Returns
    -------
    dict
        {"ocean": artist, "land": artist, "borders": artist}, any of which
        may be None if that feature is turned off. Keep this and pass it back
        in via remove_map_features() before calling this again, or the old
        artists are left behind (drawn over, not replaced).
    """

    map_template = canvas_instance.plot_characteristics_templates["map"]
    resolution = get_land_polygon_resolution(map_template["map_resolution"])
    artists = {"ocean": None, "land": None, "borders": None}

    # land/ocean colours come from the active preset unless set explicitly
    land_colour, ocean_colour = get_map_colours(map_template)

    # ocean first, so land and borders draw on top of it
    ocean_characteristics = map_template.get("ocean_polygon", {})
    if ocean_characteristics.get("visible", True):
        ocean_kwargs = {
            k: v for k, v in ocean_characteristics.items() if k != "visible"
        }
        ocean_kwargs["facecolor"] = ocean_colour
        artists["ocean"] = ax.add_feature(
            cfeature.NaturalEarthFeature(
                category="physical", name="ocean", scale=resolution, **ocean_kwargs
            )
        )

    land_kwargs = dict(map_template["land_polygon"])
    land_kwargs["facecolor"] = land_colour
    artists["land"] = ax.add_feature(
        cfeature.NaturalEarthFeature(
            category="physical",
            name="land",
            scale=resolution,
            **land_kwargs,
        )
    )

    borders_characteristics = map_template.get("borders", {})
    if borders_characteristics.get("visible", False):
        borders_kwargs = {
            k: v for k, v in borders_characteristics.items() if k != "visible"
        }
        # with no colour set, borders track the basemap's brightness - a single
        # fixed colour disappeared against the darker land/ocean options
        if not borders_kwargs.get("edgecolor"):
            borders_kwargs["edgecolor"] = map_feature_ink(map_template)
        artists["borders"] = ax.add_feature(
            cfeature.NaturalEarthFeature(
                category="cultural",
                name="admin_0_boundary_lines_land",
                scale=resolution,
                facecolor="none",
                **borders_kwargs,
            )
        )

    return artists


def remove_map_features(feature_artists):
    """
    Remove the feature artists previously returned by draw_map_features(),
    ready for it to be called again with updated settings.

    Parameters
    ----------
    feature_artists : dict
        Dict previously returned by draw_map_features().
    """

    for artist in feature_artists.values():
        if artist is not None:
            artist.remove()


def draw_map_gridlines(canvas_instance, ax, gridlines_characteristics):
    """
    Add the map's gridlines to an axis - split out from format_axis() so it can
    be re-run on its own (remove the previous Gridliner, call this again) when
    the user toggles gridlines on/off from the map settings menu, same idea as
    draw_map_features().

    Unlike land/ocean/borders, gridlines apply regardless of which map
    background is active (providentia/shaded_relief/custom image), so this
    stays a standalone function rather than folding into draw_map_features().

    Parameters
    ----------
    canvas_instance : object
        Instance of class Canvas.
    ax : cartopy.mpl.geoaxes.GeoAxes
        Map axis to draw the gridlines onto.
    gridlines_characteristics : dict
        Gridlines plot characteristics for the plot type being drawn. Passed
        in rather than read off canvas_instance, as only the dashboard keys
        its plot characteristics by base plot type - report and library pass
        the characteristics for the one plot type being made.

    Returns
    -------
    cartopy.mpl.gridliner.Gridliner or None
        The created Gridliner, or None if gridlines are turned off. Keep
        this and call .remove() on it (if not None) before calling this
        again, or the old gridlines are left behind.
    """

    if not gridlines_characteristics.get("visible", True):
        return None

    kwargs = {k: v for k, v in gridlines_characteristics.items() if k != "visible"}
    # with no colour set, gridlines track the basemap's brightness, as borders do
    if not kwargs.get("color"):
        kwargs["color"] = map_feature_ink(
            canvas_instance.plot_characteristics_templates["map"], feature="gridlines"
        )
    return ax.gridlines(crs=canvas_instance.datacrs, **kwargs)


def format_axis(
    read_instance,
    canvas_instance,
    ax,
    base_plot_type,
    plot_characteristics,
    col_ii=0,
    last_valid_row=True,
    last_row_on_page=True,
    map_extent=False,
    relevant_temporal_resolutions=None,
):
    """
    Format a plotting axis.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    canvas_instance : object
        Instance of class Canvas or Report.
    relevant_axs : list, object or dict
        Axes to format.
    relevant_data_labels : list
        Data labels corresponding to each axis.
    networkspeci : str
        Current networkspeci.
    base_plot_type : str
        Base plot type without statistical information.
    plot_type : str
        Specific plot type.
    plot_options : list
        List of plot options to apply.
    map_extent : list or bool, optional
        Map extent bounds [lonmin, lonmax, latmin, latmax] for geographic plots.
    chunk_stat : str, optional
        Chunk statistic, if relevant.
    chunk_resolution : str, optional
        Chunk resolution, if relevant.
    """

    # define default argument mutables
    if relevant_temporal_resolutions is None:
        relevant_temporal_resolutions = ["hour", "dayofweek", "month"]

    # get plot characteristics vars
    plot_characteristics_vars = list(plot_characteristics.keys())

    # get appropriate axes for nested axes
    axs_to_format = []
    temporal_resolutions_per_ax = []
    if isinstance(ax, dict):
        for relevant_temporal_resolution, sub_ax in ax.items():
            if relevant_temporal_resolution in relevant_temporal_resolutions:
                axs_to_format.append(sub_ax)
                temporal_resolutions_per_ax.append(relevant_temporal_resolution)
    elif isinstance(ax, list):
        axs_to_format = ax
        temporal_resolutions_per_ax = [""] * len(ax)
    else:
        axs_to_format.append(ax)
        temporal_resolutions_per_ax.append("")

    # iterate though relevant axes (and relevant temporal resolutions for periodic plots)
    for ax_to_format, relevant_temporal_resolution in zip(
        axs_to_format, temporal_resolutions_per_ax
    ):
        # set axis ticks and gridlines below all artists
        ax_to_format.set_axisbelow(True)

        # make axis ylabel (only on leftmost column of visible axes)?
        # if 'axis_title' in plot_characteristics_vars:
        #    ax_to_format.set_title(**plot_characteristics['axis_title'])

        # make axis xlabel?
        # if 'xlabel' in plot_characteristics_vars:
        #    ax_to_format.set_xlabel(**plot_characteristics['xlabel'])

        # make axis ylabel (only on leftmost column of visible axes)?
        # if 'ylabel' in plot_characteristics_vars:
        #    ax_to_format.set_ylabel(**plot_characteristics['ylabel'])

        # set xtick params ?
        if "xtick_params" in plot_characteristics_vars:
            ax_to_format.xaxis.set_tick_params(**plot_characteristics["xtick_params"])

        # set ytick params ?
        if "ytick_params" in plot_characteristics_vars:
            ax_to_format.yaxis.set_tick_params(**plot_characteristics["ytick_params"])

        # if are sharing xticks, and not on last row on page/last
        # valid row, then ensure current axis xticks are hidden
        if (
            ("xtick_share" in plot_characteristics_vars)
            and (not last_valid_row)
            and (not last_row_on_page)
        ):
            plt.setp(ax_to_format.get_xticklabels(), visible=False)

        # if are sharing yticks, and not on left column, then ensure current axis yticks are hidden
        if ("ytick_share" in plot_characteristics_vars) and (col_ii != 0):
            plt.setp(ax_to_format.get_yticklabels(), visible=False)

        # set xlim?
        if "xlim" in plot_characteristics_vars:
            if isinstance(plot_characteristics["xlim"], dict):
                ax_to_format.set_xlim(**plot_characteristics["xlim"])
            else:
                ax_to_format.set_xlim(plot_characteristics["xlim"])

        # set ylim?
        if "ylim" in plot_characteristics_vars:
            if isinstance(plot_characteristics["ylim"], dict):
                ax_to_format.set_ylim(**plot_characteristics["ylim"])
            else:
                ax_to_format.set_ylim(plot_characteristics["ylim"])

        # add gridlines (x and y)?
        if "grid" in plot_characteristics_vars:
            ax_to_format.grid(**plot_characteristics["grid"])

        # add x gridlines?
        if "xgrid" in plot_characteristics_vars:
            ax_to_format.xaxis.grid(**plot_characteristics["xgrid"])

        # add y gridlines?
        if "ygrid" in plot_characteristics_vars:
            ax_to_format.yaxis.grid(**plot_characteristics["ygrid"])

        # set x axis decimal places?
        if "round_decimal_places" in plot_characteristics_vars:
            if "x" in plot_characteristics["round_decimal_places"]:
                ax_to_format.xaxis.set_major_formatter(
                    ticker.FormatStrFormatter(
                        "%.{}f".format(
                            plot_characteristics["round_decimal_places"]["x"]
                        )
                    )
                )

        # set y axis decimal places?
        if "round_decimal_places" in plot_characteristics_vars:
            if "y" in plot_characteristics["round_decimal_places"]:
                ax_to_format.yaxis.set_major_formatter(
                    ticker.FormatStrFormatter(
                        "%.{}f".format(
                            plot_characteristics["round_decimal_places"]["y"]
                        )
                    )
                )

        # a density axis runs to very small numbers where what is plotted
        # covers a wide range (a station statistic like NData, or a species
        # reported in small units), and labels like "0.00005" grow wide
        # enough to push the axis label into the plot beside it, or off the
        # page - shown with one shared exponent above the axis instead
        if (
            (base_plot_type in ["distribution", "histogram"])
            and ("y" not in plot_characteristics.get("round_decimal_places", {}))
            and (ax_to_format.get_yscale() == "linear")
            and isinstance(
                ax_to_format.yaxis.get_major_formatter(), ticker.ScalarFormatter
            )
        ):
            ax_to_format.ticklabel_format(
                axis="y", style="sci", scilimits=(-3, 4), useMathText=True
            )

        # remove spines?
        if "remove_spines" in plot_characteristics_vars:
            for side in plot_characteristics["remove_spines"]:
                ax_to_format.spines[side].set_visible(False)

            for side in list(
                set(["top", "bottom", "right", "left"]).symmetric_difference(
                    plot_characteristics["remove_spines"]
                )
            ):
                ax_to_format.spines[side].set_visible(True)

        # handle formatting specific to plot types
        if base_plot_type in ["periodic", "periodic-violin"]:
            # add axis resolution label
            ax_to_format.annotate(
                canvas_instance.periodic_labels[relevant_temporal_resolution],
                **plot_characteristics["label"],
            )

            # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
            if relevant_temporal_resolution == "hour":
                plot_characteristics["xticks"] = canvas_instance.periodic_xticks[
                    relevant_temporal_resolution
                ][::3]
                ax_to_format.set_xticks(plot_characteristics["xticks"])
            else:
                plot_characteristics["xticks"] = canvas_instance.periodic_xticks[
                    relevant_temporal_resolution
                ]
                ax_to_format.set_xticks(plot_characteristics["xticks"])
                ax_to_format.set_xticklabels(
                    [
                        canvas_instance.temporal_axis_mapping_dict["short"][
                            relevant_temporal_resolution
                        ][xtick]
                        for xtick in canvas_instance.periodic_xticks[
                            relevant_temporal_resolution
                        ]
                    ]
                )

        elif base_plot_type == "map":
            # set map background

            # providentia default background
            if plot_characteristics["background"] == "providentia":
                canvas_instance.map_feature_artists = draw_map_features(
                    canvas_instance, ax_to_format
                )

            # shaded relief (cartopy default)
            elif plot_characteristics["background"] == "shaded_relief":
                ax_to_format.stock_img()

            # other type of map background
            else:
                # check file for background exists
                background_fname = join(
                    CURRENT_PATH,
                    "resources/{}.png".format(plot_characteristics["background"]),
                )
                if os.path.isfile(background_fname):
                    img = plt.imread(background_fname)
                    img_extent = (-180, 180, -90, 90)
                    ax_to_format.imshow(
                        img,
                        origin="upper",
                        extent=img_extent,
                        transform=canvas_instance.datacrs,
                    )
                else:
                    msg = "Specified map background file cannot be found."
                    show_message(read_instance, msg)

            # add gridlines ?
            if "gridlines" in plot_characteristics_vars:
                canvas_instance.map_gridliner = draw_map_gridlines(
                    canvas_instance, ax_to_format, plot_characteristics["gridlines"]
                )

            # set map extent (if wanted)
            if map_extent:
                set_map_extent(canvas_instance, ax_to_format, map_extent)

        elif base_plot_type == "fairmode-target":
            # the fixed xticks/yticks list was applied once at plot creation
            # and never adapted, so zooming into a sub-range containing none
            # of those positions left no visible ticks. Locators recompute
            # positions on every redraw, so unlike periodic this needs no
            # per-navigation hook. nbins=4 keeps the original 5-tick look
            decimal_places = plot_characteristics.get("round_decimal_places", {})
            ax_to_format.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
            ax_to_format.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
            # the diagram's x axis is signed CRMSE but conventionally
            # labelled by magnitude only - matches the original static
            # xticks, whose "labels" were already each tick's absolute
            # value (see plot_characteristics.yaml)
            ax_to_format.xaxis.set_major_formatter(
                ticker.FuncFormatter(
                    lambda value, _pos: format_tick_label(
                        abs(value), decimal_places.get("x", 2)
                    )
                )
            )
            ax_to_format.yaxis.set_major_formatter(
                ticker.FuncFormatter(
                    lambda value, _pos: format_tick_label(
                        value, decimal_places.get("y", 2)
                    )
                )
            )


def format_tick_label(value, max_decimal_places):
    """
    Format a tick value with up to a fixed number of decimal places,
    trimming any that aren't needed - e.g. 2.0 -> "2", 0.5 -> "0.5",
    0.08 -> "0.08" (at max_decimal_places=2) - rather than a fixed
    ``f"{value:.{n}f}"`` always showing every one of them (2.00, 0.50,
    ...), which looks cluttered on an axis mixing whole numbers and finer
    ones (e.g. after zooming in on part of it).

    Parameters
    ----------
    value : float
        The tick value to format.
    max_decimal_places : int
        The most decimal places to show - fewer are used if the value
        doesn't need them.

    Returns
    -------
    str
        The formatted label.
    """

    label = f"{value:.{max_decimal_places}f}"
    if "." in label:
        label = label.rstrip("0").rstrip(".")
    # rstrip can turn "-0.00" into just "-" (or "0.00" into "") - both
    # mean zero
    if label in ("", "-"):
        label = "0"

    return label


def get_no_margin_lim(ax, lim):
    """
    Get true limits of a plot area without including axis margins.

    Parameters
    ----------
    ax : object
        Matplotlib axis object to get limits from.
    lim : str
        'xlim' or 'ylim' specifying which axis to compute.

    Returns
    -------
    lower_lim : float
        Lower limit of the axis without margin.
    upper_lim : float
        Upper limit of the axis without margin.
    """

    # xlim
    if lim == "xlim":
        xlim = ax.get_xlim()
        xwidth = xlim[1] - xlim[0]
        lower_lim = xlim[0] + (0.5 * ax.margins()[0]) / (0.5 + ax.margins()[0]) * xwidth
        upper_lim = xlim[1] - (0.5 * ax.margins()[0]) / (0.5 + ax.margins()[0]) * xwidth

    # ylim
    if lim == "ylim":
        ylim = ax.get_ylim()
        ywidth = ylim[1] - ylim[0]
        lower_lim = ylim[0] + (0.5 * ax.margins()[1]) / (0.5 + ax.margins()[1]) * ywidth
        upper_lim = ylim[1] - (0.5 * ax.margins()[1]) / (0.5 + ax.margins()[1]) * ywidth

    return lower_lim, upper_lim


def get_data_lims(ax, lim, plot_options):
    """
    Get x or y limits of a plot axis based on the data actually plotted.

    Parameters
    ----------
    ax : object
        Matplotlib axis object.
    lim : str
        'xlim' or 'ylim' specifying which axis to get limits for.
    plot_options : list
        Options for the plot.

    Returns
    -------
    lower_lim : float
        Minimum value along the axis.
    upper_lim : float
        Maximum value along the axis.
    """

    # get min and max values for axis
    lines = []
    for i, line in enumerate(ax.lines):
        if lim == "xlim":
            line_data = line.get_xdata()
        elif lim == "ylim":
            line_data = line.get_ydata()
        if (list(line_data) == [0, 1]) or (
            list(line_data) == [0, 0.5] or (list(line_data) == [[1.0], [1.0]])
        ):
            continue
        lines.extend(line_data)

    # if log of an axis is active, remove values <= 0 to ensure limits are obtained correctly
    if (("logx" in plot_options) and (lim == "xlim")) or (
        ("logy" in plot_options) and (lim == "ylim")
    ):
        lines = np.array(lines)
        above_zero = lines > 0
        lines = lines[above_zero]

    # get min/max across line artists
    if len(lines) == 0:
        return np.nan, np.nan
    else:
        lower_lim = np.nanmin(lines)
        upper_lim = np.nanmax(lines)
        return lower_lim, upper_lim


def log_validity(ax, log_ax):
    """
    Determine if a log scale is valid for the given axis
    (i.e., no data ≤ 0 on that axis).

    Parameters
    ----------
    ax : object
        Matplotlib axis object.
    log_ax : str
        'logx' or 'logy', specifying which axis to check.

    Returns
    -------
    bool
        True if log scale is valid, False otherwise.
    """

    if log_ax == "logx":
        lower_lim, _ = get_data_lims(ax, "xlim", ["logx"])
    elif log_ax == "logy":
        lower_lim, _ = get_data_lims(ax, "ylim", ["logy"])

    if lower_lim < 0:
        validity = False
    else:
        validity = True

    return validity
