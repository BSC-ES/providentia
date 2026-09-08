""" Functions for plot formatting """

import copy
from datetime import datetime, timedelta
import os

import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import num2date
from matplotlib.figure import Figure
from matplotlib import ticker
from matplotlib import transforms as mtransforms
import numpy as np
import pandas as pd
from PIL import Image

from providentia.auxiliar import CURRENT_PATH, join
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

# fixed, well-known install paths for the same fonts, per OS - tried
# (potentially incomplete, e.g. inside a frozen app bundle) directory scan






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


# fixed hierarchy of calendar-aligned tick steps, finest to coarsest
# (hour multiples are every divisor of 24, so every sub-daily grid puts a
# tick on midnight; day multiples are fine-grained to give the search a
# close density match; semimonth is the 1st and 15th, which reads better
# than an arbitrary day stride over a few months)
_TIMESERIES_TICK_STEPS = [
    ("hour", 1),
    ("hour", 2),
    ("hour", 3),
    ("hour", 4),
    ("hour", 6),
    ("hour", 8),
    ("hour", 12),
    ("hour", 24),
    ("day", 1), ("day", 2), ("day", 3), ("day", 4), ("day", 5),
    ("day", 6), ("day", 7), ("day", 8), ("day", 9), ("day", 10),
    ("day", 12), ("day", 14), ("day", 16), ("day", 18), ("day", 20),
    ("day", 24), ("day", 28), ("day", 32), ("day", 36), ("day", 40),
    ("day", 45), ("day", 50), ("day", 60), ("day", 70), ("day", 80),
    ("day", 90), ("day", 100), ("day", 120), ("day", 140), ("day", 160),
    ("day", 180), ("day", 210), ("day", 240), ("day", 270), ("day", 300),
    ("day", 330), ("day", 365),
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


def _drops_a_day_start(candidates, kept_dates):
    """
    Determine if any midnight among the candidates falling strictly inside the
    span shown was left out of the kept ticks. Midnights at or outside the two
    end ticks do not count, as those are the view's own boundaries.

    Parameters
    ----------
    candidates : list
        Candidate tick datetimes
    kept_dates : list
        Tick datetimes that survived decluttering

    Returns
    -------
    bool
        True if a midnight inside the span was dropped
    """
    if not kept_dates:
        return False
    first, last = kept_dates[0], kept_dates[-1]
    kept = set(kept_dates)
    return any(
        first < candidate < last and candidate not in kept
        for candidate in candidates
        if (candidate.hour, candidate.minute, candidate.second, candidate.microsecond)
        == (0, 0, 0, 0)
    )


def _aligned_timeseries_ticks(left, right, kind, multiple):
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
        # day-ordinal aligned (day 1 is 0001-01-01), so a multi-day stride
        # always lands on the same fixed set of calendar days regardless of
        # where the visible range happens to start
        start_date = left.date()
        if (left.hour, left.minute, left.second, left.microsecond) != (
            0,
            0,
            0,
            0,
        ):
            start_date += timedelta(days=1)
        ordinal = start_date.toordinal()
        remainder = ordinal % multiple
        if remainder:
            ordinal += multiple - remainder
        step = timedelta(days=multiple)
        current = datetime.fromordinal(ordinal)
        while current <= right:
            ticks.append(current)
            current += step

    elif kind == "semimonth":
        # the 1st and 15th of every month in range - the only kind here whose
        # gaps aren't equal in real time. `multiple` is unused, but kept for a
        # consistent per-kind signature
        year, month = left.year, left.month
        while True:
            for day in (1, 15):
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
    if kind in ("day", "semimonth"):
        return dt.strftime("%Y-%m-%d")
    return dt.strftime("%Y-%m-%d %Hh")


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


def _nudge_edge_labels_onscreen(ax, renderer, min_gap_pixels):
    """
    Shift the rightmost tick label back inside the axes when it hangs off the
    visible plot area, by the minimum needed to stop its real (unpadded) box
    being clipped. The left edge is never nudged. The shift is capped at how
    far it can go before crowding the previous label, as the declutter pass
    measured every gap with this label still centred. Only meaningful for the
    final label set shown, not during the search over candidate tiers.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to nudge the labels on
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure the labels against
    min_gap_pixels : float
        Minimum gap to keep between neighbouring labels
    """

    labels = ax.xaxis.get_majorticklabels()
    n = len(labels)
    if n < 2:
        return
    ax_bbox = ax.get_window_extent(renderer)

    # every label's transform has already been reset to its un-nudged state
    # by the preceding _set_timeseries_tick_alignment() call, so only the
    # label actually being nudged here needs touching
    label = labels[-1]
    base_transform = getattr(label, "_ptv_base_transform", None)
    if base_transform is None:
        base_transform = label.get_transform()
        label._ptv_base_transform = base_transform

    # the real, unpadded box - only shift as far as is needed to stop the
    # text clipping, and not even quite that far, as a label pulled fully
    # flush with the boundary reads worse than a few px of overflow
    bbox = label.get_window_extent(renderer)
    x0 = bbox.x0
    shift_px = _edge_nudge_shift(bbox, ax_bbox)

    if shift_px and n >= 2:
        neighbour_bbox = labels[-2].get_window_extent(renderer)
        neighbour_pad = (neighbour_bbox.width * 0.15) / 2
        neighbour_x1 = neighbour_bbox.x1 + neighbour_pad
        # shift_px < 0 (moving left) - stop short of the previous
        # label's space
        min_shift = (neighbour_x1 + min_gap_pixels) - x0
        shift_px = max(shift_px, min(0, min_shift))

    if shift_px:
        shift_in = shift_px / ax.figure.dpi
        label.set_transform(
            base_transform
            + mtransforms.ScaledTranslation(shift_in, 0, ax.figure.dpi_scale_trans)
        )


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


def _disambiguate_edge_labels(left, right, kind):
    """
    Get text for the two forced start/end labels. Unlike interior ticks, the
    edges are wherever the view was zoomed or snapped to, and at a coarse
    resolution can format identically - so if the two texts come out equal,
    precision is escalated (minutes, then seconds) until they differ. Pixel
    level clash avoidance still happens in the shared declutter pass.

    Parameters
    ----------
    left : datetime.datetime
        Left edge time
    right : datetime.datetime
        Right edge time
    kind : str
        Step kind of the chosen tier

    Returns
    -------
    tuple of str
        Left and right edge label text
    """

    left_text = _edge_tick_text(left, kind)
    right_text = _edge_tick_text(right, kind)
    if left_text != right_text:
        return left_text, right_text

    if kind == "hour":
        escalations = ["%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"]
    else:
        escalations = ["%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"]
    for fmt in escalations:
        left_text, right_text = left.strftime(fmt), right.strftime(fmt)
        if left_text != right_text:
            return left_text, right_text

    # last-resort universal fallback - unreachable in practice, since
    # left < right always differ at full precision, but kept as a
    # defensive floor rather than ever returning a duplicate pair
    return (
        left.strftime("%Y-%m-%d %H:%M:%S"),
        right.strftime("%Y-%m-%d %H:%M:%S"),
    )


# fraction of the correction needed to bring a clipped right-hand edge label
# back on screen that is actually applied - a label pulled exactly flush with
# the boundary reads as more obviously shunted
_EDGE_NUDGE_FRACTION = 0.85


def _edge_nudge_shift(bbox, ax_bbox):
    """
    Get the pixel shift the rightmost label will really be given once drawn.
    The single source of truth for that number, as the space has to be reserved
    while deciding which interior ticks fit and then applied when the chosen set
    is displayed - if the two disagree, fit decisions set aside room the label
    never uses and interior ticks near the right-hand end are dropped for it.

    Parameters
    ----------
    bbox : matplotlib.transforms.Bbox
        The label's real, unpadded bounding box
    ax_bbox : matplotlib.transforms.Bbox
        The axis's bounding box

    Returns
    -------
    float
        Shift in pixels, negative, or zero when the label already fits
    """
    new_x0, _new_x1 = _onscreen_bbox_x(bbox.x0, bbox.x1, ax_bbox, "right")
    return (new_x0 - bbox.x0) * _EDGE_NUDGE_FRACTION


def _onscreen_bbox_x(x0, x1, ax_bbox, side):
    """
    Get a label's real (x0, x1) shifted just enough to stay within the axis's
    bounding box, applying the same correction _nudge_edge_labels_onscreen()
    makes once a label set is displayed - so it is accounted for as reserved
    space up front rather than discovered as a clash after the choice is made.

    Parameters
    ----------
    x0 : float
        Left edge of the label box, in pixels
    x1 : float
        Right edge of the label box, in pixels
    ax_bbox : matplotlib.transforms.Bbox
        The axis's bounding box
    side : str
        Which edge of the box might be overflowing ("left" or "right")

    Returns
    -------
    tuple of float
        Corrected (x0, x1)
    """
    if side == "left" and x0 < ax_bbox.x0:
        shift = ax_bbox.x0 - x0
        x0, x1 = x0 + shift, x1 + shift
    elif side == "right" and x1 > ax_bbox.x1:
        shift = x1 - ax_bbox.x1
        x0, x1 = x0 - shift, x1 - shift
    return x0, x1


def _measure_and_declutter(
    ax, renderer, candidates, min_gap_pixels, xlim, protected=frozenset()
):
    """
    Install the candidates as real ticks on the axis and greedily keep as many
    as fit left to right without their rendered bounding boxes crowding each
    other, always keeping the first and last. Protected candidates are taken in
    a pass of their own first, so a landmark tick is not crowded out by an
    ordinary one just before it.

    The renderer is reused as-is rather than redrawing the figure per candidate
    set, as get_window_extent() lays the string out against it fresh each call.
    xlim is re-applied straight after set_ticks(), as matplotlib silently
    expands the axis's data limits to fit any tick outside the current view,
    even with autoscale off, which would corrupt every later measurement.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to install the ticks on
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure the labels against
    candidates : list
        (datetime, text) tuples, first and last being the two edges
    min_gap_pixels : float
        Minimum gap to keep between neighbouring labels
    xlim : tuple
        Numeric x-axis limits to restore after setting ticks
    protected : frozenset, optional
        Candidate datetimes given first claim on the available space

    Returns
    -------
    tuple of list
        Kept tick datetimes and their label text, in order
    """

    candidate_dates = [c[0] for c in candidates]
    candidate_texts = [c[1] for c in candidates]
    ax.xaxis.set_ticks(candidate_dates, labels=candidate_texts)
    ax.set_xlim(*xlim)
    tick_labels = ax.xaxis.get_majorticklabels()
    _set_timeseries_tick_alignment(tick_labels)
    ax_bbox = ax.get_window_extent(renderer)

    n = len(candidates)
    label_edges = []
    for i, ((dt, _text), label) in enumerate(zip(candidates, tick_labels)):
        bbox = label.get_window_extent(renderer)
        # a small safety margin around the measured box, rather than
        # trusting it to the last pixel - two labels sitting exactly
        # min_gap_pixels apart with nothing to spare would still read
        # as touching
        pad = (bbox.width * 0.15) / 2
        x0, x1 = bbox.x0 - pad, bbox.x1 + pad
        # the right edge label is centred on a tick at the data boundary, so
        # its box can extend past the axis. Reserve exactly the shift
        # _nudge_edge_labels_onscreen() will give it, so the chosen
        # candidates already leave room. The left edge is never nudged
        if i == n - 1:
            shift = _edge_nudge_shift(bbox, ax_bbox)
            x0, x1 = x0 + shift, x1 + shift
        label_edges.append((dt, x0, x1))

    # always keep the first (left) candidate, then add later ones only if
    # they don't crowd the last kept label. `protected` candidates are taken
    # in a pass of their own first, and a plain candidate is only kept if it
    # also leaves room for the next protected one - a single greedy pass
    # would let an ordinary tick take the space a midnight needed
    interior = label_edges[1:-1]
    end_dt, end_x0, end_x1 = label_edges[-1]

    chosen_protected = []
    if protected:
        cursor_x1 = label_edges[0][2]
        for dt, x0, x1 in interior:
            if dt in protected and x0 - cursor_x1 >= min_gap_pixels:
                chosen_protected.append((dt, x0, x1))
                cursor_x1 = x1
        # the right edge is mandatory too, so give up the protected ones
        # closest to it rather than let them crowd it out
        while chosen_protected and end_x0 - chosen_protected[-1][2] < min_gap_pixels:
            chosen_protected.pop()

    kept = [label_edges[0]]
    next_protected = 0
    for dt, x0, x1 in interior:
        if next_protected < len(chosen_protected) and dt == chosen_protected[next_protected][0]:
            kept.append((dt, x0, x1))
            next_protected += 1
            continue
        if dt in protected:
            # protected, but already ruled out above - never fill its
            # place with a neighbour it would have displaced
            continue
        if x0 - kept[-1][2] < min_gap_pixels:
            continue
        if next_protected < len(chosen_protected):
            upcoming = chosen_protected[next_protected]
            if upcoming[1] - x1 < min_gap_pixels:
                continue
        kept.append((dt, x0, x1))

    # always keep the last (right) candidate too - dropping back through
    # whatever was already kept if it would otherwise crowd this one
    while len(kept) > 1:
        if end_x0 - kept[-1][2] >= min_gap_pixels:
            break
        kept.pop()
    if kept[-1][0] != end_dt:
        kept.append((end_dt, end_x0, end_x1))

    text_by_date = {dt: text for dt, text in candidates}
    kept_dates = [dt for dt, _x0, _x1 in kept]
    kept_texts = [text_by_date[dt] for dt in kept_dates]
    return kept_dates, kept_texts


def _edge_pair_fits(
    ax, renderer, edge_left, edge_right, left_text, right_text, min_gap_pixels, xlim
):
    """
    Determine whether the two edge labels fit side by side on their own,
    accounting for the same on-screen correction the right edge gets once
    displayed, so a pair is not judged as fitting only for the nudge to then
    have nowhere to go.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to install the ticks on
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure the labels against
    edge_left : datetime.datetime
        Left edge time
    edge_right : datetime.datetime
        Right edge time
    left_text : str
        Left edge label text
    right_text : str
        Right edge label text
    min_gap_pixels : float
        Minimum gap to keep between the two labels
    xlim : tuple
        Numeric x-axis limits to restore after setting ticks

    Returns
    -------
    bool
        True if the two edge labels fit
    """
    ax.xaxis.set_ticks([edge_left, edge_right], labels=[left_text, right_text])
    ax.set_xlim(*xlim)
    labels = ax.xaxis.get_majorticklabels()
    _set_timeseries_tick_alignment(labels)
    ax_bbox = ax.get_window_extent(renderer)
    left_box, right_box = (label.get_window_extent(renderer) for label in labels)
    left_pad = (left_box.width * 0.15) / 2
    right_pad = (right_box.width * 0.15) / 2
    left_x0, left_x1 = left_box.x0 - left_pad, left_box.x1 + left_pad
    right_shift = _edge_nudge_shift(right_box, ax_bbox)
    right_x0 = right_box.x0 - right_pad + right_shift
    return (right_x0 - left_x1) >= min_gap_pixels


def _shorten_edge_text(left, right, kind):
    """
    Get a shorter fallback pair of edge texts, used only when the tier's normal
    edge format does not fit the two edges side by side even with no interior
    ticks. "hour" drops the date if both edges share a day, "day"/"semimonth"
    drop the year if both share one, and "month"/"year" have nothing shorter
    that would still be meaningful.

    Parameters
    ----------
    left : datetime.datetime
        Left edge time
    right : datetime.datetime
        Right edge time
    kind : str
        Step kind of the chosen tier

    Returns
    -------
    tuple of str or None
        Shortened left and right edge text, or None if this kind has no
        shorter form to offer
    """
    if kind == "hour":
        formats = ["%Hh", "%H:%M"] if left.date() == right.date() else ["%m-%d %Hh"]
    elif kind in ("day", "semimonth") and left.year == right.year:
        formats = ["%m-%d"]
    else:
        formats = []
    for fmt in formats:
        left_text, right_text = left.strftime(fmt), right.strftime(fmt)
        if left_text != right_text:
            return left_text, right_text
    return None


def _resolve_timeseries_edges(left, right, kind, edge_aligned, data_start, data_end):
    """
    Get the datetimes to actually use for the left/right edge ticks, which are
    not always left/right themselves. A side showing the full loaded data range
    uses that true boundary, so an un-zoomed view of a full year reads
    "2018-01"/"2019-01" rather than whatever the margin padding lands on. A side
    that has been zoomed snaps inward to the nearest tick in edge_aligned, which
    is day precision (or hour precision for a sub-daily view) regardless of how
    coarse the interior ticks need to be. Never rounds outward past the true
    edge, as a tick outside the view silently widens the axis.

    Parameters
    ----------
    left : datetime.datetime
        Left view boundary
    right : datetime.datetime
        Right view boundary
    kind : str
        Step kind of the chosen tier
    edge_aligned : list
        Day or hour precision aligned ticks to snap a zoomed edge to
    data_start : datetime.datetime or None
        True start of the loaded data, None if not known
    data_end : datetime.datetime or None
        True end of the loaded data, None if not known

    Returns
    -------
    tuple of datetime.datetime
        Left and right edge times
    """
    if data_start is not None and left <= data_start:
        edge_left = data_start
    elif edge_aligned:
        edge_left = edge_aligned[0]
    else:
        edge_left = left

    if data_end is not None and right >= data_end:
        edge_right = data_end
    elif edge_aligned:
        edge_right = edge_aligned[-1]
    else:
        edge_right = right

    # collision guard: a single-element aligned list (or a data range
    # narrower than one tick step) could pick the same point for both
    if edge_left >= edge_right:
        edge_left, edge_right = left, right

    return edge_left, edge_right


def _fit_timeseries_ticks(
    ax, renderer, left, right, kind, aligned, min_gap_pixels, data_start, data_end, xlim
):
    """
    Build the label candidates for one tier (the edges plus whatever aligned
    ticks fall between them) and find how many actually fit, reporting both how
    many interior candidates were offered and how many were kept so the caller
    can tell whether this tier's density suits the available space.

    A left-to-right greedy keep is used only to find out how many interior ticks
    the space can hold, as it can leave a gap-toothed result. A second attempt
    then picks that many at an even stride across the full offered set, and is
    used whenever it fits at least as many, so what is shown reads as an
    intentional coarser resolution rather than an arbitrary subset of a finer one.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to install the ticks on
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure the labels against
    left : datetime.datetime
        Left view boundary
    right : datetime.datetime
        Right view boundary
    kind : str
        Step kind of this tier
    aligned : list
        Aligned tick datetimes offered by this tier
    min_gap_pixels : float
        Minimum gap to keep between neighbouring labels
    data_start : datetime.datetime or None
        True start of the loaded data, None if not known
    data_end : datetime.datetime or None
        True end of the loaded data, None if not known
    xlim : tuple
        Numeric x-axis limits to restore after setting ticks

    Returns
    -------
    tuple
        Kept tick datetimes, their label text, the number of interior
        candidates offered and the number kept
    """

    # edges always snap at day precision, or whole-hour precision for a
    # sub-daily view - and to every hour, not this tier's own multiple, so an
    # edge sits as close to where the view really begins and ends as a
    # sensible unit allows, rather than being dragged onto the interior grid
    edge_aligned = _aligned_timeseries_ticks(left, right, "hour", 1) if kind == "hour" else _aligned_timeseries_ticks(left, right, "day", 1)
    edge_left, edge_right = _resolve_timeseries_edges(
        left, right, kind, edge_aligned, data_start, data_end
    )
    left_text, right_text = _disambiguate_edge_labels(edge_left, edge_right, kind)
    if not _edge_pair_fits(
        ax, renderer, edge_left, edge_right, left_text, right_text, min_gap_pixels, xlim
    ):
        # the two edges, both always shown, don't fit side by side at their
        # normal per-tier precision - shorten just enough for the pair to
        # fit, keeping the original text if this kind has nothing shorter
        shortened = _shorten_edge_text(edge_left, edge_right, kind)
        if shortened is not None:
            left_text, right_text = shortened
    interior = [
        (dt, _format_timeseries_tick(dt, kind))
        for dt in aligned
        if edge_left < dt < edge_right
    ]
    # drop an interior candidate showing identical text to the edge next to
    # it (e.g. a mid-month right edge and a month-start tick both reading
    # "2018-05") - the edge wins, as it is the one guaranteed to stay
    if interior and interior[0][1] == left_text:
        interior = interior[1:]
    if interior and interior[-1][1] == right_text:
        interior = interior[:-1]

    n_offered = len(interior)
    edges = [(edge_left, left_text), (edge_right, right_text)]
    if n_offered == 0:
        kept_dates, kept_texts = _measure_and_declutter(
            ax, renderer, edges, min_gap_pixels, xlim
        )
        return kept_dates, kept_texts, 0, 0

    # on a sub-daily view the start of a day is the only landmark there is,
    # so midnights get first claim on the space. Above day resolution every
    # candidate is already a day boundary, so there is nothing to single out
    protected = (
        frozenset(
            dt
            for dt, _text in interior
            if (dt.hour, dt.minute, dt.second, dt.microsecond) == (0, 0, 0, 0)
        )
        if kind == "hour"
        else frozenset()
    )

    all_candidates = [edges[0]] + interior + [edges[1]]
    kept_dates, kept_texts = _measure_and_declutter(
        ax, renderer, all_candidates, min_gap_pixels, xlim, protected
    )
    n_kept = len(kept_dates) - 2
    if n_kept == n_offered:
        return kept_dates, kept_texts, n_offered, n_kept  # everything fit

    if n_kept <= 0:
        return kept_dates, kept_texts, n_offered, 0

    if protected:
        # the protected pass above already placed the day boundaries and
        # filled around them; re-picking an evenly-spaced subset below
        # would choose purely by spacing again and undo exactly that
        return kept_dates, kept_texts, n_offered, n_kept

    # an even-stride subset reads as an intentional coarser resolution rather
    # than a leftover from the greedy pass. Greedy's own count is a lower
    # bound, not the maximum, so every size from n_offered down to 1 is a
    # candidate - fitting is monotonic in size, so the largest clean size is
    # found with a binary search rather than by trying every one
    def _even_stride_attempt(size):
        # stride across the full edge-to-edge span, not just the interior
        # list's own span - striding the interior alone always anchors the
        # first and last chosen point immediately next to an edge, however
        # sparse the size
        idx = sorted(
            {
                min(n_offered - 1, max(0, round(k * (n_offered + 1) / (size + 1)) - 1))
                for k in range(1, size + 1)
            }
        )
        chosen_interior = [interior[i] for i in idx]
        candidates = [edges[0]] + chosen_interior + [edges[1]]
        dates, texts = _measure_and_declutter(
            ax, renderer, candidates, min_gap_pixels, xlim
        )
        return dates, texts, len(dates) - 2

    best_even = None  # (dates, texts, size) - largest size confirmed clean
    lo, hi = 1, n_offered
    while lo <= hi:
        mid = (lo + hi) // 2
        dates, texts, kept_at_mid = _even_stride_attempt(mid)
        if kept_at_mid == mid:
            best_even = (dates, texts, mid)
            lo = mid + 1
        else:
            hi = mid - 1

    if best_even is not None and best_even[2] >= n_kept:
        even_dates, even_texts, even_n_kept = best_even
        return even_dates, even_texts, n_offered, even_n_kept
    return kept_dates, kept_texts, n_offered, n_kept


def compute_timeseries_xticks(
    ax, left, right, max_ticks=6, min_gap_pixels=12, data_start=None, data_end=None,
    data_resolution_seconds=None,
):
    """
    Set "nice" x-axis tick positions and labels for a timeseries date
    range directly on `ax`, always labelling both edges - with ticks in
    between snapped to a fixed hierarchy of calendar-aligned resolutions
    (see _TIMESERIES_TICK_STEPS) rather than a generic locator's full
    range of possible steps, so every tick lands on a boundary that
    actually means something (never e.g. a half hour).

    The edges themselves are not always `left`/`right` exactly - see
    _resolve_timeseries_edges(): a side still showing the full loaded
    data range reads as that range's own clean start/end (e.g. the
    default view of a full year reads "2018-01"/"2019-01", not a
    margin-padded value a few percent past it), while a side that's
    actually been zoomed in snaps to the nearest sensible tick instead
    of an exact, often visually clunky boundary.

    Density adapts to the axis's actual, current pixel width rather
    than a fixed guess at how many ticks "should" fit: tiers are tried
    from finest to coarsest, and for each candidate tier that isn't
    already ruled out on count alone, its labels are really installed,
    drawn, and measured (_fit_timeseries_ticks) to see whether they all
    survive decluttering untouched - the first (finest, most
    informative) tier where nothing had to be dropped is what's used, so
    a wide panel naturally ends up with more, closer-together ticks than
    a narrow one showing the same span, without either ever clashing.

    Three earlier versions of this decluttering step measured
    candidates *before* they were actually on screen - a fraction of
    the total data range, then a prediction of each label's rendered
    pixel width via get_text_width_height_descent() - and still let
    labels clash on a real machine despite passing the same kind of
    check in a test environment: a predicted width is only as good as
    the font/DPI assumptions behind it, and those can differ (font
    substitution, HiDPI scaling, hinting) between wherever this gets
    tested and where it actually runs. So this measures real,
    already-rendered label bounding boxes instead - whatever
    font/DPI/renderer is in play on the machine actually running this
    is exactly what gets measured, by construction.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The real axis to set ticks on - candidates are installed on it
        and it is actually drawn so each label's real rendered bounding
        box can be read back.
    left : datetime.datetime
        Start of the visible x-axis range.
    right : datetime.datetime
        End of the visible x-axis range.
    max_ticks : int, default 6
        Ceiling on how many *interior* ticks a tier may offer to even be
        considered - the two edges are always additionally shown on top
        of this, so the real total can run up to max_ticks + 2. Not a
        target: the tier actually used is whichever fits this ceiling
        *and* the axis's real available width (see above), so the
        result is commonly fewer than max_ticks on a narrow panel and
        can be noticeably more on a wide one.
    min_gap_pixels : float, default 12
        Minimum gap, in pixels, required between two adjacent labels'
        edges for both to be kept.
    data_start, data_end : datetime.datetime, optional
        The true, configured start/end of the loaded data (not the
        current view) - see _resolve_timeseries_edges(). Left as None,
        both edges are always treated as zoomed (the same behaviour as
        before this parameter existed).
    data_resolution_seconds : float, optional
        The loaded data's own sampling interval, in seconds (e.g. 3600
        for hourly). When the current view is zoomed to no wider than
        this, there is at most one real observation actually in it, so
        a single tick at that observation's own time is shown instead of
        two edge labels either side of a gap with nothing real in it.
        Left as None, this collapse never happens (the same behaviour as
        before this parameter existed).

    Returns
    -------
    xticks : list of datetime.datetime
        The tick positions actually set on `ax`, sorted, always
        starting and ending at whatever _resolve_timeseries_edges()
        settled on for that side (unless the two are so close together
        that even their own labels can't help but overlap - both are
        still kept, since having a start/end label at all wins over
        avoiding that one unavoidable clash).
    """

    full_precision = "%Y-%m-%d %H:%M:%S"

    if left >= right:
        dates = [left, right] if left < right else [left]
        ax.xaxis.set_ticks(dates, labels=[d.strftime(full_precision) for d in dates])
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
            return [snapped]

    # one real draw, to settle the axis's own layout - reused as the
    # renderer for every candidate set tried below instead of redrawing
    # the whole figure (data, other axes, everything) each time; see
    # _measure_and_declutter() for why that's still an accurate measure
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

    # try tiers finest to coarsest, skipping any whose interior tick count
    # exceeds max_ticks, then score each surviving tier on how many ticks it
    # kept, with a 2x handicap for calendar-anchored kinds - a month start is
    # a landmark a reader recognises, where a large day stride lands on a
    # date with no significance. Highest score anywhere wins

    # judge the count gate below against where the edges will actually end
    # up, not the raw view bounds - a side showing the full loaded range
    # resolves to data_start/data_end, which is narrower than left/right (an
    # autoscaled view sits a margin outside the real data). Counting against
    # the padded bounds counts that margin's own extra month as an interior
    # candidate, which tipped e.g. a Jan-to-July view's "month,1" over
    # max_ticks it would otherwise have cleared
    count_left = data_start if (data_start is not None and left <= data_start) else left
    count_right = data_end if (data_end is not None and right >= data_end) else right

    best = None
    best_score = -1
    best_calendar_kind = False
    for step_kind, multiple in _TIMESERIES_TICK_STEPS:
        if step_kind != "hour" and best is not None and best[3] == 0:
            # every "hour" tier kept zero interior candidates, and coarser
            # steps can only have less to offer, so stop rather than measure
            # all ~50 remaining day/month/year steps
            break
        step_ticks = _aligned_timeseries_ticks(left, right, step_kind, multiple)
        interior_count = sum(1 for t in step_ticks if count_left < t < count_right)
        if interior_count > max_ticks:
            continue
        result = _fit_timeseries_ticks(
            ax, renderer, left, right, step_kind, step_ticks, min_gap_pixels,
            data_start, data_end, xlim,
        )
        _kept_dates, _kept_texts, n_offered, n_kept = result
        # a coarser tier with nothing to offer trivially keeps 0 too, but
        # must never outscore a finer tier that had something to show -
        # guaranteed, as a 0-kept tier scores 0 whatever the handicap
        calendar_kind = step_kind in ("semimonth", "month", "year")
        score = n_kept * (2 if calendar_kind else 1)
        if step_kind == "hour" and _drops_a_day_start(step_ticks, _kept_dates):
            # midnight is always among an "hour" tier's candidates, but the
            # thinning picks by spacing alone and readily drops it. Losing a
            # day start is what makes a sub-daily view read as random, so
            # score it far below any midnight-preserving option - still
            # proportional to n_kept, so the fullest wins when none can
            score *= 0.01
        # a tie goes to the calendar-anchored tier - day tiers are tried
        # first, so a plain ">" would let a same-scoring day tier win purely
        # for being found first, undoing the calendar handicap
        if best is None or score > best_score or (
            score == best_score
            and score > 0
            and calendar_kind
            and not best_calendar_kind
        ):
            best = result
            best_score = score
            best_calendar_kind = calendar_kind

    if best is None:
        # every tier offered more interior candidates than max_ticks allows
        # - stride evenly through the coarsest rather than showing an
        # unbounded number of ticks
        step_kind, multiple = _TIMESERIES_TICK_STEPS[-1]
        step_ticks = _aligned_timeseries_ticks(left, right, step_kind, multiple)
        stride = max(1, len(step_ticks) // max(1, max_ticks))
        best = _fit_timeseries_ticks(
            ax, renderer, left, right, step_kind, step_ticks[::stride], min_gap_pixels,
            data_start, data_end, xlim,
        )

    kept_dates, kept_texts, _n_offered, _n_kept = best
    ax.xaxis.set_ticks(kept_dates, labels=kept_texts)
    ax.set_xlim(*xlim)
    _set_timeseries_tick_alignment(ax.xaxis.get_majorticklabels())
    _nudge_edge_labels_onscreen(ax, renderer, min_gap_pixels)

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
                # always label the visible start and end, with ticks between
                # snapped to a hierarchy of calendar-aligned resolutions (see
                # compute_timeseries_xticks()) instead of evenly slicing the
                # range into a fixed number of pieces. data_start/data_end
                # let a side still showing the full range read as that
                # range's own boundary rather than a margin-padded value,
                # each side judged independently
                max_ticks = plot_characteristics["xtick_alteration"]["max_ticks"]
                data_start = _parse_yyyymmdd(getattr(read_instance, "start_date", None))
                data_end = _parse_yyyymmdd(getattr(read_instance, "end_date", None))
                active_resolution = getattr(read_instance, "active_resolution", None) or getattr(
                    read_instance, "resolution", None
                )
                data_resolution_seconds = _RESOLUTION_SECONDS.get(active_resolution)
                for ax in relevant_axs_active:
                    compute_timeseries_xticks(
                        ax, left, right, max_ticks=max_ticks,
                        data_start=data_start, data_end=data_end,
                        data_resolution_seconds=data_resolution_seconds,
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
                    canvas_instance, ax_to_format
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
