""" Functions and classes that allow interactivity of dashboard plots """

import copy
import datetime

from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.widgets import _SelectorWidget
import numpy as np
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QApplication, QToolTip, QWidget

from .dashboard_elements import (
    set_formatting,
    set_cursor,
    unset_cursor,
    LegendInlineEditor,
    LegendEditorCommitFilter,
)
from .plot_aux import get_display_label, get_map_extent, get_hex_code


class LassoSelector(_SelectorWidget):
    """
    Selection curve of an arbitrary shape.
    For the selector to remain responsive you must keep a reference to it.
    The selected path can be used in conjunction with `~.Path.contains_point`
    to select data points from an image.
    In contrast to `Lasso`, `LassoSelector` is written with an interface
    similar to `RectangleSelector` and `SpanSelector`, and will continue to
    interact with the Axes until disconnected.
    Example usage::
        ax = plt.subplot()
        ax.plot(x, y)
        def onselect(verts):
            print(verts)
        lasso = LassoSelector(ax, onselect)
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The parent Axes for the widget.
    onselect : function
        Whenever the lasso is released, the *onselect* function is called and
        passed the vertices of the selected path.
    useblit : bool, default: True
        Whether to use blitting for faster drawing (if supported by the
        backend). See the tutorial :doc:`/tutorials/advanced/blitting`
        for details.
    props : dict, optional
        Properties with which the line is drawn, see `matplotlib.lines.Line2D`
        for valid properties. Default values are defined in ``mpl.rcParams``.
    button : `.MouseButton` or list of `.MouseButton`, optional
        The mouse buttons used for rectangle selection.  Default is ``None``,
        which corresponds to all buttons.
    """

    def __init__(self, ax, onselect, useblit=True, props=None, button=None):
        """
        Initialise class

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes
        onselect : callable
            Callback function called when a selection is completed.
            It receives the list of polygon vertices.
        useblit : bool
            Whether to use blitting for faster drawing
        props : dict
            Line2D properties for the selection outline (e.g., color, linewidth).
        button : int or list of int
            Mouse button(s) used to start the selection.
        """

        super().__init__(ax, onselect, useblit=useblit, button=button)
        self.verts = None
        props = {
            **(props if props is not None else {}),
            # Note that self.useblit may be != useblit, if the canvas doesn't
            # support blitting.
            "animated": self.useblit,
            "visible": False,
        }
        line = Line2D([], [], **props)
        self.ax.add_line(line)
        self._selection_artist = line

        return None

    def _press(self, event):
        """
        Handle mouse press event to start a new selection

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse event containing the click position
        """

        self.verts = [self._get_data(event)]
        self._selection_artist.set_visible(True)

        return None

    def _onmove(self, event):
        """
        Handle mouse movement during selection

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse move event used to update the selection path
        """

        if self.verts is None:
            return
        self.verts.append(self._get_data(event))
        self._selection_artist.set_data(list(zip(*self.verts)))

        self.update()

        return None

    def _release(self, event):
        """
        Handle mouse release event and finalise the selection

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse release event marking the end of the selection
        """

        if self.verts is not None:
            self.verts.append(self._get_data(event))
            self.onselect(self.verts)
        self._selection_artist.set_data([[], []])
        self._selection_artist.set_visible(False)
        self.verts = None

        return None


# opacity the legend editor draws its text at - Qt rasterises the same font
# at the same size about 15% heavier than matplotlib's Agg, and heavier again
# on macOS. Lower this if the text reads heavier than the label
_EDITOR_INK_MATCH = 0.78

# Qt font families registered from matplotlib's own font files, keyed by
# file path - see _qt_font_matching(). Registering the same file twice is
# harmless but pointless, and the lookup happens on every legend rename.
_REGISTERED_FONT_FAMILIES = {}


def _qt_font_matching(font_properties):
    """
    A QFont backed by the very font file matplotlib resolved for these
    properties, so the two lay text out the same way.

    Naming a family and hoping Qt picks the same face isn't enough: Qt
    resolves family names through its own substitution table and readily
    lands on a different face with the same name, or a fallback with
    similar overall proportions but different individual glyph widths.
    Calibrating the total string width papers over that, but the caret
    sits at a *character* offset - so the per-glyph disagreement still
    accumulates left to right, which is the drift this fixes. Loading
    matplotlib's actual file removes the disagreement at source.

    Parameters
    ----------
    font_properties : matplotlib.font_manager.FontProperties
        The properties matplotlib is rendering the text with.

    Returns
    -------
    QtGui.QFont
        Font using the same face, at its default size - the caller sets
        the size.
    """

    try:
        font_path = font_manager.findfont(font_properties)
    except Exception:
        font_path = None

    family = _REGISTERED_FONT_FAMILIES.get(font_path) if font_path else None
    if font_path and family is None:
        try:
            font_id = QtGui.QFontDatabase.addApplicationFont(font_path)
            families = QtGui.QFontDatabase.applicationFontFamilies(font_id)
            family = families[0] if families else None
        except Exception:
            family = None
        _REGISTERED_FONT_FAMILIES[font_path] = family

    if family is None:
        # nothing loadable - fall back to the declared family name, which
        # is the best guess available and what this did before
        declared = font_properties.get_family()
        family = declared[0] if declared else ""

    return QtGui.QFont(family)


def _calibrate_editor_font(
    editor_font,
    text,
    renderer,
    font_properties,
    dpi_ratio,
    nominal_size,
    min_size,
    max_size,
    logical_dpi_y,
):
    """
    Set `editor_font`'s size and letter spacing so Qt's caret lands
    where matplotlib actually drew each character, and return the pixel
    size chosen.

    Two things have to be right. The obvious one is the size, fitted here
    by matching Qt's width for this exact string to matplotlib's, with
    whatever is left over spread across the characters as letter spacing.

    The subtler one is hinting. By default Qt grid-fits glyphs to whole
    pixels when it renders, so each character's *drawn* advance is rounded
    even though QFontMetricsF reports the true fractional value - and the
    caret is placed using the rounded ones. Fractions of a pixel per
    character accumulate into a caret that is visibly right of where it
    should be by the end of a long label, while every measurement taken
    of it still says the two agree to within a pixel. matplotlib does no
    such rounding, so the editor is switched to unhinted outline
    rendering to match it.

    Parameters
    ----------
    editor_font : QtGui.QFont
        Font to configure, already using matplotlib's own font file (see
        _qt_font_matching()).
    text : str
        The label being edited, used as the calibration sample.
    renderer : matplotlib.backend_bases.RendererBase
        Renderer to measure matplotlib's own text with.
    font_properties : matplotlib.font_manager.FontProperties
        The properties matplotlib draws the label with.
    dpi_ratio : float
        The canvas' device pixel ratio - matplotlib measures in device
        pixels, the widget is laid out in logical ones.
    nominal_size : int
        Starting estimate for the pixel size.
    min_size, max_size : int
        Bounds the size is kept within, in logical pixels.
    logical_dpi_y : float
        The widget's logical vertical DPI, for converting the wanted
        pixel height into the fractional point size Qt takes.

    Returns
    -------
    float
        The size set on the font, in logical pixels.
    """

    # unhinted, outline-rendered: keeps Qt's drawn advances fractional,
    # the way matplotlib's are
    editor_font.setHintingPreference(QtGui.QFont.PreferNoHinting)
    # grayscale antialiasing, as matplotlib's Agg uses. Left alone Qt takes
    # the platform's subpixel smoothing, which on macOS thickens stems
    editor_font.setStyleStrategy(QtGui.QFont.NoSubpixelAntialias)

    # sized in fractional points, not whole pixels - setPixelSize only takes
    # an integer, so matching the width meant rounding up and pulling the
    # extra back out with letter spacing, leaving the glyphs drawn too large
    logical_dpi = max(logical_dpi_y, 1.0)

    def set_size(logical_pixels):
        editor_font.setPointSizeF(max(logical_pixels, 1.0) * 72.0 / logical_dpi)

    target_pixels = min(max(float(nominal_size), min_size), max_size)
    set_size(target_pixels)
    if not text:
        return target_pixels

    try:
        target_width = (
            renderer.get_text_width_height_descent(text, font_properties, False)[0]
            / max(dpi_ratio, 0.1)
        )
    except Exception:
        return target_pixels
    if target_width <= 0:
        return target_pixels

    qt_width = QtGui.QFontMetricsF(editor_font).horizontalAdvance(text)
    if qt_width > 0:
        # self-correcting: absorbs any mismatch in the dpi/pixel-ratio
        # assumption behind nominal_size, and bounded so a bad measurement
        # can't run away into an editor big enough to cover the window
        target_pixels = min(
            max(target_pixels * target_width / qt_width, min_size), max_size
        )
        set_size(target_pixels)

    # only the sub-pixel remainder is left to letter spacing now, rather
    # than a whole pixel of over-size per character
    qt_width = QtGui.QFontMetricsF(editor_font).horizontalAdvance(text)
    if qt_width > 0:
        spacing = (target_width - qt_width) / len(text)
        limit = target_pixels * 0.4
        editor_font.setLetterSpacing(
            QtGui.QFont.AbsoluteSpacing, min(max(spacing, -limit), limit)
        )

    return target_pixels


def _editor_ink_offset(editor):
    """
    Where `editor` actually puts the first pixel of its text, relative to
    its own top-left corner.

    A QLineEdit does not draw its text at its origin: it keeps an internal
    horizontal margin (measured here at 3px, and style- and
    platform-dependent), and centres the line vertically by font metrics
    rather than by ink. Positioning the widget by its corner therefore put
    the text a few pixels right of, and slightly below, the label it was
    standing in for. Rather than assume any of those constants, the widget
    is rendered once offscreen and its ink located directly, so the caller
    can line that ink up with matplotlib's.

    Parameters
    ----------
    editor : QtWidgets.QLineEdit
        The editor, already carrying its final font, text and size.

    Returns
    -------
    tuple of (float, float) or None
        (x, y) of the first inked pixel within the widget, or None if
        nothing could be measured.
    """

    try:
        image = QtGui.QImage(editor.size(), QtGui.QImage.Format_ARGB32)
        image.fill(QtCore.Qt.transparent)
        editor.render(image)
        buffer = image.bits()
        buffer.setsize(image.byteCount())
        pixels = np.frombuffer(buffer, np.uint8).reshape(
            image.height(), image.width(), 4
        )
        inked = pixels[..., 3] > 40
        if not inked.any():
            return None
        rows = np.where(inked.any(axis=1))[0]
        columns = np.where(inked.any(axis=0))[0]
        return float(columns.min()), float(rows.min())
    except Exception:
        return None


def zoom_map_func(canvas_instance, event):
    """
    Handle scroll-wheel zooming on the map axis

    Parameters
    ----------
    canvas_instance : object
        Canvas instance
    event : matplotlib.backend_bases.MouseEvent
        Scroll event providing mouse position, scroll direction, and axis info
    """

    if event.inaxes == canvas_instance.plot_axes["map"]:
        # lock canvas drawing if can, else return
        if canvas_instance.figure.canvas.widgetlock.locked():
            if not canvas_instance.figure.canvas.widgetlock.isowner(
                canvas_instance.zoom_map_event
            ):
                return None
        else:
            canvas_instance.figure.canvas.widgetlock(canvas_instance.zoom_map_event)

        # get the current x and y limits
        current_xlim = canvas_instance.plot_axes["map"].get_xlim()
        current_ylim = canvas_instance.plot_axes["map"].get_ylim()

        # get position of cursor
        xdata = event.xdata
        ydata = event.ydata
        base_scale = canvas_instance.plot_characteristics["map"]["base_scale"]

        if event.button == "up":
            # deal with zoom in
            scale_factor = base_scale
        elif event.button == "down":
            # deal with zoom out
            scale_factor = 1 / base_scale
        else:
            # exceptions
            scale_factor = 1

        if event.button == "up" or event.button == "down":
            # set new limits
            canvas_instance.plot_axes["map"].set_xlim(
                [
                    xdata - (xdata - current_xlim[0]) / scale_factor,
                    xdata + (current_xlim[1] - xdata) / scale_factor,
                ]
            )
            canvas_instance.plot_axes["map"].set_ylim(
                [
                    ydata - (ydata - current_ylim[0]) / scale_factor,
                    ydata + (current_ylim[1] - ydata) / scale_factor,
                ]
            )

            # save map extent (in data coords)
            canvas_instance.read_instance.map_extent = get_map_extent(canvas_instance)

            # re-derive automatic marker size/opacity for the new zoom
            # level (a no-op if automatic sizing is off)
            canvas_instance.apply_automatic_marker_style()

            # draw changes
            canvas_instance.figure.canvas.draw_idle()

            # update buttons (previous-forward) history
            canvas_instance.read_instance.navi_toolbar.push_current()
            canvas_instance.read_instance.navi_toolbar.set_history_buttons()

        # unlock canvas drawing
        if canvas_instance.figure.canvas.widgetlock.isowner(
            canvas_instance.zoom_map_event
        ):
            canvas_instance.figure.canvas.widgetlock.release(
                canvas_instance.zoom_map_event
            )

    return None


def picker_block_func(canvas_instance, event):
    """
    Enable or disable legend picking depending on where the click occurs

    Parameters
    ----------
    canvas_instance : object
        Canvas instance
    event : matplotlib.backend_bases.MouseEvent
        Mouse event used to determine which axis was clicked
    """

    if event.inaxes == canvas_instance.plot_axes["legend"]:
        # unblock legend picker in legend
        canvas_instance.lock_legend_pick = False

    else:
        # block legend picker
        canvas_instance.lock_legend_pick = True

    return None


def legend_picker_func(canvas_instance, event):
    """
    Handle clicks on legend items to toggle visibility of plotted data

    Parameters
    ----------
    canvas_instance : object
        Canvas instance
    event : matplotlib.backend_bases.PickEvent
        Pick event providing the clicked legend artist
    """

    # pick_event fires for any pickable artist in the figure, not just the
    # legend (the map's station scatter has its own picker), so ignore
    # anything that isn't a legend text
    if not isinstance(event.artist, Text):
        return None

    # get legend label information - gid carries the real data label
    # regardless of what's actually displayed (see get_display_label() in
    # plotting.py), falling back to the displayed text itself if for any
    # reason gid isn't set
    legend_label = event.artist
    data_label = legend_label.get_gid() or legend_label.get_text()

    # double-click renames the legend/display text for this data label,
    # instead of toggling its visibility
    if event.mouseevent.dblclick:
        # matplotlib fires a normal (non-double) pick_event for the *first*
        # click of a double-click too, which would otherwise have already
        # queued a toggle below - cancel it so double-clicking to rename
        # doesn't also hide that data label's plots until clicked again
        pending_timer = getattr(canvas_instance, "_pending_legend_toggle_timer", None)
        if pending_timer is not None:
            pending_timer.stop()
            canvas_instance._pending_legend_toggle_timer = None
        rename_legend_label(canvas_instance, legend_label, data_label)
        return None

    # defer the single-click toggle by the system's double-click interval, in
    # case a second click arrives and this turns out to be the first half of a
    # double-click - stopped above if so
    timer = QtCore.QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(
        lambda: _toggle_legend_visibility(canvas_instance, legend_label, data_label)
    )
    canvas_instance._pending_legend_toggle_timer = timer
    timer.start(QApplication.instance().doubleClickInterval())

    return None


def _toggle_legend_visibility(canvas_instance, legend_label, data_label):
    """
    Show/hide a data label's plotted elements - the actual single-click
    legend behaviour, called after legend_picker_func()'s short debounce
    confirms the click wasn't the first half of a double-click.

    Parameters
    ----------
    canvas_instance : object
        Canvas instance
    legend_label : matplotlib.text.Text
        The clicked legend text artist.
    data_label : str
        The data label's real identifier (from legend_label's gid - see
        legend_picker_func()).
    """

    if not canvas_instance.lock_legend_pick:
        if canvas_instance.plot_elements:
            # lock legend pick
            canvas_instance.lock_legend_pick = True

            if data_label not in canvas_instance.plot_elements["data_labels_active"]:
                visible = True
                # put observations label always first in pop-ups on hover
                if data_label == canvas_instance.read_instance.observations_data_label:
                    canvas_instance.plot_elements["data_labels_active"].insert(
                        0, data_label
                    )
                # put model labels in the same order as in the legend
                else:
                    canvas_instance.plot_elements["data_labels_active"].insert(
                        list(canvas_instance.read_instance.experiments.values()).index(
                            data_label
                        )
                        + 1,
                        data_label,
                    )
            else:
                visible = False
                canvas_instance.plot_elements["data_labels_active"].remove(data_label)

            # iterate through plot types stored in plot_elements (if have selected stations)
            if len(canvas_instance.relative_selected_station_inds) > 0:
                for plot_type in canvas_instance.plot_elements:
                    if plot_type not in [
                        "data_labels_active",
                        "metadata",
                        "map",
                        "heatmap",
                        "table",
                        "statsummary",
                        "boxplot",
                    ]:
                        # get currently selected options for plot
                        plot_options = canvas_instance.current_plot_options[plot_type]

                        # get active (absolute / bias)
                        active = canvas_instance.plot_elements[plot_type]["active"]

                        # change visibility of plot elements (if data label in plot elements dictionary)
                        if (
                            data_label
                            in canvas_instance.plot_elements[plot_type][active]
                        ):
                            for element_type in canvas_instance.plot_elements[
                                plot_type
                            ][active][data_label]:
                                for plot_element in canvas_instance.plot_elements[
                                    plot_type
                                ][active][data_label][element_type]:
                                    if visible:
                                        plot_element.set_visible(True)
                                    else:
                                        plot_element.set_visible(False)

                # the boxplot's categories - position, spacing and tick
                # labels alike - depend on which data labels are currently
                # shown, so it is redrawn from scratch here rather than
                # having the hidden one's box merely toggled invisible in
                # place, which left its own tick and label behind (see
                # make_boxplot() in plotting.py)
                if "boxplot" in canvas_instance.read_instance.active_dashboard_plots:
                    canvas_instance.update_associated_active_dashboard_plot("boxplot")

            # change font weight of label
            legend_label._fontproperties = canvas_instance.legend.get_texts()[
                0
            ]._fontproperties.copy()
            if visible:
                legend_label.set_fontweight("bold")
            else:
                legend_label.set_fontweight("regular")

            # draw changes
            canvas_instance.figure.canvas.draw_idle()

            # unlock legend pick
            canvas_instance.lock_legend_pick = False

    return None


def rename_legend_label(canvas_instance, legend_label, data_label):
    """
    Edit a data label's display name in place, overlaying a QLineEdit on top of
    the double-clicked legend text - Enter or clicking away commits, Escape
    cancels. The new name is shown in the legend and any other plot drawing data
    labels as text (see get_display_label() in plot_aux.py), but data_label
    itself, the real identifier used for data selection and style lookups, is
    never touched. The override is cleared by emptying the field, and by the
    next data load.

    Parameters
    ----------
    canvas_instance : object
        Canvas instance
    legend_label : matplotlib.text.Text
        The double-clicked legend text artist, used to position the editor
        directly over it.
    data_label : str
        The data label's real identifier (from legend_label's gid - see
        legend_picker_func()).
    """

    read_instance = canvas_instance.read_instance
    current_display = read_instance.legend_label_overrides.get(data_label, data_label)

    # the QLineEdit captures keystrokes, cursor and selection, and draws the
    # text itself while the matplotlib label is hidden - see the note on
    # editor.setFont() below for why Qt draws both rather than matplotlib
    figure_canvas = canvas_instance.figure.canvas
    renderer = figure_canvas.get_renderer()
    bbox = legend_label.get_window_extent(renderer)
    dpi_ratio = figure_canvas.devicePixelRatioF()
    fig_height = canvas_instance.figure.bbox.height
    figure_dpi = canvas_instance.figure.dpi
    # the legend's Text is anchored left/baseline, so its position is exactly
    # the pen origin matplotlib drew from - a better reference than the ink
    # bounding box, which shifts with whichever letters the label contains
    # and put the editor slightly down and to the right of the label
    pen_x, baseline_y = legend_label.get_transform().transform(
        legend_label.get_position()
    )
    x = pen_x / dpi_ratio
    mpl_baseline_from_top = (fig_height - baseline_y) / dpi_ratio

    # match the concrete font matplotlib resolves to, then scale its point
    # size so Qt's rendered width of this exact string agrees with
    # matplotlib's - the caret is positioned from Qt's own glyph advances, so
    # a narrower font lands it short of where the visible text ends
    font_properties = legend_label.get_fontproperties()
    editor_font = _qt_font_matching(font_properties)
    editor_font.setBold(font_properties.get_weight() in ("bold", "heavy", 700, 800, 900))
    editor_font.setItalic(font_properties.get_style() in ("italic", "oblique"))

    # size the editor's font by matching its rendered width to matplotlib's
    # for this exact string. As _qt_font_matching() has Qt using the same
    # font file, the two agree on each glyph rather than merely the total,
    # which is what the caret position depends on
    mpl_width = bbox.width / dpi_ratio
    declared_size = font_properties.get_size()
    # a legend label is always ordinary UI-sized text, so anything outside
    # this range is a bad measurement - dpi and device pixel ratio can both
    # be stale right after the window moves to a differently scaled screen
    MIN_EDITOR_PIXELS, MAX_EDITOR_PIXELS = 6, 40

    nominal_size = int(round(declared_size * figure_dpi / 72.0 / max(dpi_ratio, 0.1)))
    nominal_size = min(max(nominal_size, MIN_EDITOR_PIXELS), MAX_EDITOR_PIXELS)
    pixel_size = _calibrate_editor_font(
        editor_font,
        current_display,
        renderer,
        font_properties,
        dpi_ratio,
        nominal_size,
        MIN_EDITOR_PIXELS,
        MAX_EDITOR_PIXELS,
        float(figure_canvas.logicalDpiY()),
    )

    original_text = legend_label.get_text()

    editor = LegendInlineEditor(figure_canvas)
    editor.setText(current_display)
    editor.setFont(editor_font)
    editor.setFrame(False)
    editor.setTextMargins(0, 0, 0, 0)
    # the editor draws its own text, in the legend's colour, and the
    # matplotlib label is hidden while it does. Qt drawing both the text and
    # the caret is the only way the two cannot disagree - the caret sits
    # between the letters by construction rather than by calibration. The
    # size comes from _calibrate_editor_font(), fitted to the width
    # matplotlib measures for this exact string
    editor.setFont(editor_font)

    # height comes after the match, since that is what settles the final
    # font size. Bounded by that size rather than by any measured bbox -
    # deriving it from a measurement that may itself be wrong is no
    # protection, which is how the runaway editor came back once before.
    widget_height = QtGui.QFontMetricsF(editor_font).height()
    widget_height = min(max(widget_height, MIN_EDITOR_PIXELS), pixel_size * 3)

    editor_colour = QtGui.QColor(get_hex_code(legend_label.get_color()))
    editor.setStyleSheet(
        "QLineEdit {"
        "  border: none; background: transparent; padding: 0;"
        f"  color: rgba({editor_colour.red()}, {editor_colour.green()}, "
        f"{editor_colour.blue()}, {_EDITOR_INK_MATCH});"
        "}"
    )

    def editor_width_for(text):
        """
        Get the width the editor needs to hold the given text. Qt will not put
        the caret past the right-hand edge of the widget, so a fixed width also
        caps how far into the text the caret can be placed - sizing from the
        text itself, with room for a few more characters, removes that limit.

        Parameters
        ----------
        text : str
            Text the editor has to hold

        Returns
        -------
        int
            Width in pixels
        """

        metrics = QtGui.QFontMetricsF(editor_font)
        return int(round(metrics.horizontalAdvance(text or " ") + metrics.height() * 2))

    # Qt centres the text line vertically in the widget, so its baseline
    # sits this far below the widget's top - line the two baselines up
    # rather than guessing from box centres.
    editor_metrics = QtGui.QFontMetricsF(editor_font)
    qt_baseline_from_top = (
        widget_height + editor_metrics.ascent() - editor_metrics.descent()
    ) / 2

    # keep the editor wholly inside the canvas whatever the measurements
    # said - the last line of defence against a stale bbox or pixel ratio
    # putting it somewhere absurd
    canvas_width = max(figure_canvas.width(), 1)
    canvas_height = max(figure_canvas.height(), 1)
    editor_x = int(round(min(max(x, 0), canvas_width - 1)))
    editor_y = int(round(mpl_baseline_from_top - qt_baseline_from_top))
    editor_y = int(min(max(editor_y, 0), canvas_height - 1))
    editor_width = min(editor_width_for(current_display), canvas_width - editor_x)
    editor.setGeometry(editor_x, editor_y, editor_width, int(round(widget_height)))

    # with the size settled, line the editor's own ink up horizontally with
    # the ink matplotlib drew, absorbing the widget's internal text margin
    ink_offset = _editor_ink_offset(editor)
    if ink_offset is not None:
        ink_dx, _ = ink_offset
        target_ink_x = bbox.x0 / dpi_ratio
        editor_x = int(round(min(max(target_ink_x - ink_dx, 0), canvas_width - 1)))
        editor_width = min(editor_width_for(current_display), canvas_width - editor_x)
        editor.setGeometry(editor_x, editor_y, editor_width, int(round(widget_height)))

    # vertically it is the baselines that are lined up, not the ink. A text's
    # extent in matplotlib is reported from the font's metrics rather than
    # from the glyphs the string actually holds, so its top sits at the
    # ascender line whether or not anything reaches it - matching it to ink
    # lifted a name of short letters ("cams") by the height of the ascenders
    # it does not have. The caret's own box gives Qt's line position, and is
    # the same whatever the letters, so the two agree for any name
    caret = editor.cursorRect()
    qt_baseline_from_top = caret.top() + QtGui.QFontMetricsF(editor_font).ascent()
    editor_y = int(
        round(
            min(max(mpl_baseline_from_top - qt_baseline_from_top, 0), canvas_height - 1)
        )
    )
    editor.move(editor.x(), editor_y)
    editor.show()
    editor.setFocus()
    # not selectAll() - an active selection has no blinking caret, and as the
    # selection highlight is styled transparent above, a full selection left
    # nothing visible until an arrow key collapsed it back to a caret
    editor.end(False)

    # keep a reference so PyQt doesn't garbage-collect the Python wrapper
    # out from under a still-alive (Qt-parented) widget before it's used
    canvas_instance._legend_inline_editor = editor

    # the label is hidden for the edit and the canvas redrawn once without
    # it, so the editor's text is the only rendering of the name on screen -
    # nothing to keep in step per keystroke, and no redraw while typing
    legend_label.set_visible(False)
    figure_canvas.draw()

    def sync_preview(new_text):
        # only the widget's width needs maintaining, so the caret can
        # always reach the end of what has been typed (editor_width_for())
        editor.resize(
            min(editor_width_for(new_text), canvas_width - editor_x), editor.height()
        )

    editor.textChanged.connect(sync_preview)

    # editingFinished fires on Enter, but clicking elsewhere on the canvas is
    # handled by matplotlib's own pick machinery rather than Qt
    # click-to-focus, so the editor doesn't reliably lose focus - hence the
    # commit_filter below. This guard keeps those paths from running twice
    handled = {"done": False}

    def cleanup():
        editor.textChanged.disconnect(sync_preview)
        # whichever way the edit ends, the hidden artist has to come back
        # - update_legend() rebuilds the legend on commit, but nothing
        # does on cancel
        legend_label.set_visible(True)
        QApplication.instance().removeEventFilter(commit_filter)
        editor.deleteLater()

    def commit():
        if handled["done"]:
            return
        handled["done"] = True
        cleanup()

        # Providentia logo busy cursor, same as every other action in the
        # app that takes a perceptible moment - see set_cursor()/
        # unset_cursor() in dashboard_elements.py
        read_instance.cursor_function = set_cursor(
            read_instance.cursor_function, "rename_legend_label"
        )
        # force the cursor change to paint before the work below restores it,
        # as the set+unset can otherwise happen within one event loop pass.
        # ExcludeUserInputEvents so this pump cannot process a fresh click
        QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)

        new_text = editor.text().strip()
        if (not new_text) or (new_text == data_label):
            read_instance.legend_label_overrides.pop(data_label, None)
        else:
            read_instance.legend_label_overrides[data_label] = new_text

        # update_legend() rebuilds the legend from scratch (new Text
        # objects), so the live-edited legend_label above is discarded
        # rather than needing to be reset here
        canvas_instance.update_legend()

        # a rename doesn't touch any underlying data or selection, and
        # statsummary's table and the boxplot's category labels are the only
        # other places a display label appears - the broader
        # update_associated_active_dashboard_plots() re-fetches station data
        # and redraws every active plot, which felt slow. The boxplot redraw
        # also re-decides whether its labels now fit (see
        # fit_boxplot_xticklabels()), as a new name can be shorter or longer
        # than the one it last measured
        for plot_type in ("statsummary", "boxplot"):
            if plot_type in read_instance.active_dashboard_plots:
                canvas_instance.update_associated_active_dashboard_plot(plot_type)

        # draw(), not draw_idle(): a deferred repaint lands after the
        # cursor below has been restored, leaving the slow part of the
        # rename showing the plain pointer instead of the Providentia one
        canvas_instance.figure.canvas.draw()

        unset_cursor(read_instance.cursor_function, "rename_legend_label")

    def cancel():
        if handled["done"]:
            return
        handled["done"] = True
        cleanup()
        # nothing rebuilds the legend on cancel, so put the original text
        # back rather than leaving whatever was last typed on screen
        legend_label.set_text(original_text)
        canvas_instance.figure.canvas.draw_idle()

    editor.editingFinished.connect(commit)
    editor.escapePressed.connect(cancel)

    commit_filter = LegendEditorCommitFilter(editor, commit)
    QApplication.instance().installEventFilter(commit_filter)

    return None


class HoverAnnotation(object):
    def __init__(self, canvas_instance):
        """Initialise class

        Parameters
        ----------
        canvas_instance : object
            Canvas instance
        """

        self.canvas_instance = canvas_instance

        # set up formatting for canvas annotations
        self.canvas_instance.figure.canvas = set_formatting(
            self.canvas_instance.figure.canvas,
            self.canvas_instance.read_instance.formatting_dict["canvas_annotation"],
        )
        self.canvas_instance.figure.canvas.setToolTip("")
        QToolTip.hideText()

        # set up formatting for canvas annotation vline
        self.canvas_instance.canvas_annotation_vline = set_formatting(
            QWidget(self.canvas_instance),
            self.canvas_instance.read_instance.formatting_dict[
                "canvas_annotation_vline"
            ],
        )
        self.canvas_instance.canvas_annotation_vline.hide()

        return None

    def pointer_over_menu(self, event):
        """
        Whether the pointer is over one of the settings menus, which are Qt
        widgets sat on top of the canvas rather than anything drawn into it.

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse move event

        Returns
        -------
        bool
            Whether a menu is under the pointer
        """

        canvas_instance = self.canvas_instance
        if (event.x is None) or (event.y is None):
            return False

        # matplotlib counts from the bottom left of the figure in its own
        # pixels; Qt counts from the top left of the widget in screen ones
        ratio = canvas_instance.figure.canvas.devicePixelRatioF()
        position = QtCore.QPoint(
            int(event.x / ratio),
            int((canvas_instance.figure.bbox.height - event.y) / ratio),
        )

        # the hover line and the legend's editor are the canvas's own, and
        # follow the pointer around - they are not something to keep clear of
        not_menus = (
            getattr(canvas_instance, "canvas_annotation_vline", None),
            getattr(canvas_instance, "_legend_inline_editor", None),
        )

        for child in canvas_instance.children():
            if (
                isinstance(child, QWidget)
                and (child not in not_menus)
                and child.isVisible()
                and child.geometry().contains(position)
            ):
                return True

        return False

    def hover_legend_label(self, event):
        """
        Show the pointer over a legend label as editable, and say how.

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse move event

        Returns
        -------
        bool
            Whether a legend label is being hovered over
        """

        legend = getattr(self.canvas_instance, "legend", None)
        figure_canvas = self.canvas_instance.figure.canvas

        hovered = None
        if (legend is not None) and (event.inaxes is not None):
            if event.inaxes == self.canvas_instance.plot_axes.get("legend"):
                for legend_label in legend.texts:
                    if legend_label.get_visible() and legend_label.contains(event)[0]:
                        hovered = legend_label
                        break

        if hovered is None:
            # only put the pointer back if this is what changed it, so that
            # nothing else setting a cursor is undone from here
            if getattr(self.canvas_instance, "_legend_hover_cursor", False):
                figure_canvas.unsetCursor()
                self.canvas_instance._legend_hover_cursor = False
            return False

        if not getattr(self.canvas_instance, "_legend_hover_cursor", False):
            figure_canvas.setCursor(QtCore.Qt.IBeamCursor)
            self.canvas_instance._legend_hover_cursor = True
        figure_canvas.setToolTip("Double-click to rename")

        return True

    def hover_annotation(self, event):
        """
        Handle hover events to display point annotations on timeseries, scatter, distribution, taylor
        and FAIRMODE target plots

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse move event used to determine hovered axis and point
        """

        # hide annotation and vline
        self.canvas_instance.figure.canvas.setToolTip("")
        QToolTip.hideText()
        self.canvas_instance.canvas_annotation_vline.hide()

        # an open settings menu sits over the canvas, and its controls are
        # there to be used - a tooltip following the pointer across them
        # gets in the way of reading and clicking them
        if self.pointer_over_menu(event):
            return None

        # a legend label can be renamed by double-clicking it, which nothing
        # on screen says, so hovering one shows the text cursor and a hint -
        # the same invitation an editable field on a page gives
        if self.hover_legend_label(event):
            return None

        # identify which axis is currently being hovered over
        plot_type = None
        for test_plot_type in self.canvas_instance.plot_axes:
            if not plot_type:
                if test_plot_type in ["periodic", "periodic-violin"]:
                    if hasattr(
                        self.canvas_instance.read_instance,
                        "periodic_relevant_temporal_resolutions",
                    ):
                        for (
                            resolution
                        ) in (
                            self.canvas_instance.read_instance.periodic_relevant_temporal_resolutions
                        ):
                            if (
                                event.inaxes
                                == self.canvas_instance.plot_axes[test_plot_type][
                                    resolution
                                ]
                            ):
                                plot_type = copy.deepcopy(test_plot_type)
                                break
                elif test_plot_type == "fairmode-statsummary":
                    for i in range(len(self.canvas_instance.plot_axes[test_plot_type])):
                        if (
                            event.inaxes
                            == self.canvas_instance.plot_axes[test_plot_type][i]
                        ):
                            plot_type = copy.deepcopy(test_plot_type)
                            break
                else:
                    if event.inaxes == self.canvas_instance.plot_axes[test_plot_type]:
                        plot_type = copy.deepcopy(test_plot_type)
                        break

        # if an axis is being hovered over then now check if a point is being hovered over
        if plot_type:
            # add active axis to self
            if plot_type in ["periodic", "periodic-violin"]:
                self.ax = self.canvas_instance.plot_axes[plot_type][resolution]
            elif plot_type == "fairmode-statsummary":
                self.ax = self.canvas_instance.plot_axes[plot_type][i]
            else:
                self.ax = self.canvas_instance.plot_axes[plot_type]

            # activate hover over plot
            if plot_type == "map":
                search_plot = "stations_scatter"
            elif plot_type == "periodic":
                search_plot = "periodic_plots"
            elif plot_type == "periodic-violin":
                search_plot = "violin_plot"
            else:
                search_plot = "{}_plot".format(plot_type.replace("-", "_"))

            # get plot element name
            if plot_type == "periodic":
                plot_element_name = "plot_{}".format(resolution)
            elif plot_type == "periodic-violin":
                plot_element_name = "Median_plot_{}".format(resolution)
            else:
                plot_element_name = "plot"

            if (hasattr(self.canvas_instance.plotting, search_plot)) and (
                plot_type in self.canvas_instance.plot_elements
            ):
                is_contained = False

                for data_label in self.canvas_instance.plot_elements[
                    "data_labels_active"
                ]:
                    # skip observations for bias plot
                    if (
                        (
                            plot_type
                            in [
                                "timeseries",
                                "distribution",
                                "histogram",
                                "periodic",
                                "periodic-violin",
                            ]
                        )
                        and (
                            self.canvas_instance.plot_elements[plot_type]["active"]
                            == "bias"
                        )
                        and (
                            data_label
                            == self.canvas_instance.read_instance.observations_data_label
                        )
                    ):
                        continue

                    # do not annotate if plot is cleared
                    if (
                        data_label
                        not in self.canvas_instance.plot_elements[plot_type][
                            self.canvas_instance.plot_elements[plot_type]["active"]
                        ].keys()
                    ):
                        continue

                    if plot_type == "map":
                        (
                            is_contained,
                            annotation_index,
                        ) = self.canvas_instance.plotting.stations_scatter.contains(
                            event
                        )
                    else:
                        # do no annotate if hidedata is active
                        if (
                            len(
                                self.canvas_instance.plot_elements[plot_type][
                                    self.canvas_instance.plot_elements[plot_type][
                                        "active"
                                    ]
                                ][data_label][plot_element_name]
                            )
                            == 0
                        ):
                            continue
                        line = self.canvas_instance.plot_elements[plot_type][
                            self.canvas_instance.plot_elements[plot_type]["active"]
                        ][data_label][plot_element_name]
                        for n_point, point in enumerate(line):
                            is_contained, annotation_index = point.contains(event)
                            if is_contained:
                                if plot_type in [
                                    "fairmode-target",
                                    "fairmode-statsummary",
                                ]:
                                    annotation_index = {
                                        "ind": np.array([n_point], dtype=np.int32)
                                    }
                                break

                    if is_contained:
                        break

                # point is being hovered over?
                if is_contained:
                    # add event coordinates to self (handling pixel scaling)
                    self.x = round(
                        event.x / self.canvas_instance.read_instance.devicePixelRatio()
                    )
                    self.y = round(
                        event.y / self.canvas_instance.read_instance.devicePixelRatio()
                    )

                    # add xdata event to self
                    self.xdata = event.xdata

                    # update annotation and show vline
                    func = getattr(
                        self, "update_{}_annotation".format(plot_type.replace("-", "_"))
                    )
                    if plot_type in ["periodic", "periodic-violin"]:
                        func(annotation_index, resolution)
                    elif plot_type in [
                        "fairmode-target",
                        "fairmode-statsummary",
                        "scatter",
                        "taylor",
                    ]:
                        func(annotation_index, data_label)
                    else:
                        func(annotation_index)

        return None

    def update_map_annotation(self, annotation_index):
        """
        Update the tooltip annotation for a hovered station on the map

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered station
        """

        # retrieve stations references and coordinates
        station_name = self.canvas_instance.read_instance.station_names[
            self.canvas_instance.read_instance.networkspeci
        ][
            self.canvas_instance.active_map_valid_station_inds[
                annotation_index["ind"][0]
            ]
        ]
        station_reference = self.canvas_instance.read_instance.station_references[
            self.canvas_instance.read_instance.networkspeci
        ][
            self.canvas_instance.active_map_valid_station_inds[
                annotation_index["ind"][0]
            ]
        ]
        station_location = self.canvas_instance.plotting.stations_scatter.get_offsets()[
            annotation_index["ind"][0]
        ]
        station_value = self.canvas_instance.z_statistic[annotation_index["ind"][0]]

        # create annotation text
        text_label = ("Station: {0}\n").format(station_name)
        text_label += ("Reference: {0}\n").format(station_reference)
        text_label += ("Longitude: {0:.2f}\n").format(station_location[0])
        text_label += ("Latitude: {0:.2f}\n").format(station_location[1])
        text_label += ("{0}: {1:.{2}f}").format(
            self.canvas_instance.map_z_stat.currentText(),
            station_value,
            self.canvas_instance.plot_characteristics["map"][
                "marker_annotate_rounding"
            ],
        )

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_timeseries_annotation(self, annotation_index):
        """
        Update the tooltip annotation for a hovered point on the timeseries

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        """

        # initialise annotation text
        text_label = ""

        # iterate through active data labels
        for data_label in self.canvas_instance.plot_elements["data_labels_active"]:
            # skip observations for bias plot
            if (
                self.canvas_instance.plot_elements["timeseries"]["active"] == "bias"
                and data_label
                == self.canvas_instance.read_instance.observations_data_label
            ):
                continue

            # do not annotate if plot is cleared
            if (
                data_label
                not in self.canvas_instance.plot_elements["timeseries"][
                    self.canvas_instance.plot_elements["timeseries"]["active"]
                ].keys()
            ):
                continue

            # retrieve time and concentration
            line = self.canvas_instance.plot_elements["timeseries"][
                self.canvas_instance.plot_elements["timeseries"]["active"]
            ][data_label]["plot"][0]
            time = line.get_xdata()[annotation_index["ind"][0]]
            concentration = line.get_ydata()[annotation_index["ind"][0]]

            # first valid data label?
            if not text_label:
                # update vline position
                self.update_vline_position()

                # add time to annotation text
                text_label += ("<p style='white-space:pre'><i>Time: {0}</i>").format(
                    time.astype("datetime64[us]")
                    .astype(datetime.datetime)
                    .strftime("%d/%m/%Y %H:%M:%S")
                )

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label][
                "colour"
            ]

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
                hex_colour,
                get_display_label(self.canvas_instance.read_instance, data_label),
                concentration,
                self.canvas_instance.plot_characteristics["timeseries"][
                    "marker_annotate_rounding"
                ],
            )

        # end formatting of text label
        text_label += "</p>"

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None

    def update_scatter_annotation(self, annotation_index, data_label):
        """
        Update the tooltip annotation for a hovered point on the scatter

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        data_label : str
            Data label
        """

        # initialise annotation text
        text_label = ""

        # do not annotate if plot is cleared
        if (
            data_label
            not in self.canvas_instance.plot_elements["scatter"][
                self.canvas_instance.plot_elements["scatter"]["active"]
            ].keys()
        ):
            return None

        # retrieve concentrations in x and y axis
        line = self.canvas_instance.plot_elements["scatter"][
            self.canvas_instance.plot_elements["scatter"]["active"]
        ][data_label]["plot"][0]
        x = line.get_xdata()[annotation_index["ind"][0]]
        y = line.get_ydata()[annotation_index["ind"][0]]

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label][
            "colour"
        ]

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        text_label += ('<font color="{0}">{1}</font>').format(
            hex_colour,
            get_display_label(self.canvas_instance.read_instance, data_label),
        )
        # observations label
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
            hex_colour,
            "x",
            x,
            self.canvas_instance.plot_characteristics["scatter"][
                "marker_annotate_rounding"
            ],
        )
        # model label
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
            hex_colour,
            "y",
            y,
            self.canvas_instance.plot_characteristics["scatter"][
                "marker_annotate_rounding"
            ],
        )

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_fairmode_target_annotation(self, annotation_index, data_label):
        """
        Update the tooltip annotation for a hovered point on the FAIRMODE target plot

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        data_label : str
            Data label
        """

        # initialise annotation text
        text_label = ""

        # do not annotate if plot is cleared
        if (
            data_label
            not in self.canvas_instance.plot_elements["fairmode-target"][
                self.canvas_instance.plot_elements["fairmode-target"]["active"]
            ].keys()
        ):
            return None

        # retrieve CRMSE / β·RMSᵤ and Mean Bias / β·RMSᵤ
        line = self.canvas_instance.plot_elements["fairmode-target"][
            self.canvas_instance.plot_elements["fairmode-target"]["active"]
        ][data_label]["plot"][annotation_index["ind"][0]]
        x = line.get_xdata()[0]
        y = line.get_ydata()[0]

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label][
            "colour"
        ]

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        text_label += ('<font color="{0}">{1}</font>').format(
            hex_colour,
            get_display_label(self.canvas_instance.read_instance, data_label),
        )
        # CRMSE
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
            hex_colour,
            "CRMSE / β·RMSᵤ",
            x,
            self.canvas_instance.plot_characteristics["fairmode-target"][
                "marker_annotate_rounding"
            ],
        )
        # MB
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
            hex_colour,
            "MB / β·RMSᵤ",
            y,
            self.canvas_instance.plot_characteristics["fairmode-target"][
                "marker_annotate_rounding"
            ],
        )

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_fairmode_statsummary_annotation(self, annotation_index, data_label):
        """
        Update the tooltip annotation for a hovered point on the FAIRMODE statistics summary plot

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        data_label : str
            Data label
        """

        # initialise annotation text
        text_label = ""

        # do not annotate if plot is cleared
        if (
            data_label
            not in self.canvas_instance.plot_elements["fairmode-statsummary"][
                self.canvas_instance.plot_elements["fairmode-statsummary"]["active"]
            ].keys()
        ):
            return None

        # retrieve CRMSE / β·RMSᵤ and Mean Bias / β·RMSᵤ
        line = self.canvas_instance.plot_elements["fairmode-statsummary"][
            self.canvas_instance.plot_elements["fairmode-statsummary"]["active"]
        ][data_label]["plot"][annotation_index["ind"][0]]

        # get closest value to the point currently hovered
        x = line.get_xdata()
        closest_value = min(x, key=lambda i: abs(i - self.xdata))

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label][
            "colour"
        ]

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        display_label = get_display_label(
            self.canvas_instance.read_instance, data_label
        )
        text_label += f'<font color="{hex_colour}">{display_label}</font>'

        # CRMSE
        rounding = self.canvas_instance.plot_characteristics["fairmode-statsummary"][
            "marker_annotate_rounding"
        ]
        text_label += (
            f'<br><font color="{hex_colour}">x: {closest_value:.{rounding}f}</font>'
        )

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_distribution_annotation(self, annotation_index):
        """
        Update the tooltip annotation for a hovered point on the distribution plot

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        """

        # initialise annotation text
        text_label = ""

        # iterate through active data labels
        for data_label in self.canvas_instance.plot_elements["data_labels_active"]:
            # skip observations for bias plot
            if (
                self.canvas_instance.plot_elements["distribution"]["active"] == "bias"
            ) and (
                data_label == self.canvas_instance.read_instance.observations_data_label
            ):
                continue

            # do not annotate if plot is cleared
            if (
                data_label
                not in self.canvas_instance.plot_elements["distribution"][
                    self.canvas_instance.plot_elements["distribution"]["active"]
                ].keys()
            ):
                continue

            # retrieve concentration and density
            line = self.canvas_instance.plot_elements["distribution"][
                self.canvas_instance.plot_elements["distribution"]["active"]
            ][data_label]["plot"][0]
            concentration = line.get_xdata()[annotation_index["ind"][0]]
            density = line.get_ydata()[annotation_index["ind"][0]]

            # first valid data label?
            if not text_label:
                # update vline position
                self.update_vline_position()

                # create annotation text
                text_label += (
                    "<p style='white-space:pre'><i>{0}: {1:.{2}f}</i>"
                ).format(
                    self.canvas_instance.read_instance.species[0],
                    concentration,
                    self.canvas_instance.plot_characteristics["distribution"][
                        "marker_annotate_rounding"
                    ],
                )

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label][
                "colour"
            ]

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
                hex_colour,
                get_display_label(self.canvas_instance.read_instance, data_label),
                density,
                self.canvas_instance.plot_characteristics["distribution"][
                    "marker_annotate_rounding"
                ],
            )

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None

    def update_histogram_annotation(self, annotation_index):
        """
        Update the tooltip annotation for a hovered bin on the histogram plot

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        """

        # initialise annotation text
        text_label = ""

        # iterate through active data labels
        for data_label in self.canvas_instance.plot_elements["data_labels_active"]:
            # skip observations for bias plot
            if (
                self.canvas_instance.plot_elements["histogram"]["active"] == "bias"
            ) and (
                data_label == self.canvas_instance.read_instance.observations_data_label
            ):
                continue

            # do not annotate if plot is cleared
            if (
                data_label
                not in self.canvas_instance.plot_elements["histogram"][
                    self.canvas_instance.plot_elements["histogram"]["active"]
                ].keys()
            ):
                continue

            # retrieve the bin's edge and density
            line = self.canvas_instance.plot_elements["histogram"][
                self.canvas_instance.plot_elements["histogram"]["active"]
            ][data_label]["plot"][0]
            concentration = line.get_xdata()[annotation_index["ind"][0]]
            density = line.get_ydata()[annotation_index["ind"][0]]

            # first valid data label?
            if not text_label:
                # update vline position
                self.update_vline_position()

                # create annotation text
                text_label += (
                    "<p style='white-space:pre'><i>{0}: {1:.{2}f}</i>"
                ).format(
                    self.canvas_instance.read_instance.species[0],
                    concentration,
                    self.canvas_instance.plot_characteristics["histogram"][
                        "marker_annotate_rounding"
                    ],
                )

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label][
                "colour"
            ]

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
                hex_colour,
                get_display_label(self.canvas_instance.read_instance, data_label),
                density,
                self.canvas_instance.plot_characteristics["histogram"][
                    "marker_annotate_rounding"
                ],
            )

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None

    def update_taylor_annotation(self, annotation_index, data_label):
        """
        Update the tooltip annotation for a hovered point on the Taylor diagram plot

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        """

        # initialise annotation text
        text_label = ""

        # do not annotate if plot is cleared
        if (
            data_label
            not in self.canvas_instance.plot_elements["taylor"][
                self.canvas_instance.plot_elements["taylor"]["active"]
            ].keys()
        ):
            return None

        # retrieve time and concentration
        line = self.canvas_instance.plot_elements["taylor"][
            self.canvas_instance.plot_elements["taylor"]["active"]
        ][data_label]["plot"][0]
        corr_stat = line.get_xdata()[annotation_index["ind"][0]]
        stddev = line.get_ydata()[annotation_index["ind"][0]]

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label][
            "colour"
        ]

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        text_label += ('<font color="{0}">{1}</font>').format(
            hex_colour,
            get_display_label(self.canvas_instance.read_instance, data_label),
        )
        # corr stat
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
            hex_colour,
            self.canvas_instance.plot_characteristics["taylor"]["corr_stat"],
            np.cos(corr_stat),
            self.canvas_instance.plot_characteristics["taylor"][
                "marker_annotate_rounding"
            ],
        )
        # stddev
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
            hex_colour,
            "StdDev",
            stddev,
            self.canvas_instance.plot_characteristics["taylor"][
                "marker_annotate_rounding"
            ],
        )

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_periodic_annotation(self, annotation_index, resolution):
        """
        Update the tooltip annotation for a hovered point on the periodic plot

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        resolution : str
            Temporal resolution
        """

        # initialise annotation text
        text_label = ""

        # iterate through active data labels
        for data_label in self.canvas_instance.plot_elements["data_labels_active"]:
            # skip observations for bias plot
            if (
                self.canvas_instance.plot_elements["periodic"]["active"] == "bias"
                and data_label
                == self.canvas_instance.read_instance.observations_data_label
            ):
                continue

            # do not annotate if plot is cleared
            if (
                data_label
                not in self.canvas_instance.plot_elements["periodic"][
                    self.canvas_instance.plot_elements["periodic"]["active"]
                ].keys()
            ):
                continue

            # retrieve time and concentration
            line = self.canvas_instance.plot_elements["periodic"][
                self.canvas_instance.plot_elements["periodic"]["active"]
            ][data_label]["plot_" + resolution][0]
            time = line.get_xdata()[annotation_index["ind"][0]]
            concentration = line.get_ydata()[annotation_index["ind"][0]]

            # first valid data label?
            if not text_label:
                # update vline position
                self.update_vline_position()

                # create annotation text
                if resolution == "hour":
                    resolution_text = "Hour"
                    time_text = time
                else:
                    time_options = [
                        self.canvas_instance.temporal_axis_mapping_dict["long"][
                            resolution
                        ][xtick]
                        for xtick in self.canvas_instance.periodic_xticks[resolution]
                    ]
                    if resolution == "dayofweek":
                        time_text = time_options[time]
                        resolution_text = "Day"
                    elif resolution == "month":
                        time_text = time_options[time - 1]
                        resolution_text = "Month"
                text_label += ("<p style='white-space:pre'><i>{0}: {1}</i>").format(
                    resolution_text, time_text
                )

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label][
                "colour"
            ]

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
                hex_colour,
                get_display_label(self.canvas_instance.read_instance, data_label),
                concentration,
                self.canvas_instance.plot_characteristics["periodic"][
                    "marker_annotate_rounding"
                ],
            )

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None

    def update_periodic_violin_annotation(self, annotation_index, resolution):
        """
        Update the tooltip annotation for a hovered point on the periodic violin plot

        Parameters
        ----------
        annotation_index : dict
            Dictionary containing the index of the hovered point
        resolution : str
            Temporal resolution
        """

        # initialise annotation text
        text_label = ""

        # iterate through active data labels
        for data_label in self.canvas_instance.plot_elements["data_labels_active"]:
            # skip observations for bias plot
            if (
                self.canvas_instance.plot_elements["periodic-violin"]["active"]
                == "bias"
                and data_label
                == self.canvas_instance.read_instance.observations_data_label
            ):
                continue

            # do not annotate if plot is cleared
            if (
                data_label
                not in self.canvas_instance.plot_elements["periodic-violin"][
                    self.canvas_instance.plot_elements["periodic-violin"]["active"]
                ].keys()
            ):
                continue

            # retrieve time and concentration
            line = self.canvas_instance.plot_elements["periodic-violin"][
                self.canvas_instance.plot_elements["periodic-violin"]["active"]
            ][data_label]["Median_plot_" + resolution][0]
            time = line.get_xdata()[annotation_index["ind"][0]]
            concentration = line.get_ydata()[annotation_index["ind"][0]]

            # first valid data label?
            if not text_label:
                # update vline position
                self.update_vline_position()

                # create annotation text
                if resolution == "hour":
                    resolution_text = "Hour"
                    time_text = time
                else:
                    time_options = [
                        self.canvas_instance.temporal_axis_mapping_dict["long"][
                            resolution
                        ][xtick]
                        for xtick in self.canvas_instance.periodic_xticks[resolution]
                    ]
                    if resolution == "dayofweek":
                        time_text = time_options[time]
                        resolution_text = "Day"
                    elif resolution == "month":
                        time_text = time_options[time - 1]
                        resolution_text = "Month"
                text_label += ("<p style='white-space:pre'><i>{0}: {1}</i>").format(
                    resolution_text, time_text
                )

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label][
                "colour"
            ]

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(
                hex_colour,
                get_display_label(self.canvas_instance.read_instance, data_label),
                concentration,
                self.canvas_instance.plot_characteristics["periodic-violin"][
                    "marker_annotate_rounding"
                ],
            )

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None

    def update_vline_position(self):
        """
        Update the position of the vertical line on hover
        """

        # get current canvas width / height
        (
            canvas_width,
            canvas_height,
        ) = self.canvas_instance.figure.canvas.get_width_height()

        # get axis ylim (in data coordinates)
        ylim = self.ax.get_ylim()

        # transform matplotlib ylim data coordinates to display coordinates, handling pixel scaling
        ymin_display_mpl = round(
            self.ax.transData.transform([0, ylim[0]])[1]
            / self.canvas_instance.read_instance.devicePixelRatio()
        )
        ymax_display_mpl = round(
            self.ax.transData.transform([0, ylim[1]])[1]
            / self.canvas_instance.read_instance.devicePixelRatio()
        )

        # transform matplotlib display coordinates to Qt display coordinates
        ymin_display_qt = int(canvas_height - ymin_display_mpl)
        ymax_display_qt = int(canvas_height - ymax_display_mpl)

        # calculate length of vline to plot
        vline_length = ymin_display_qt - ymax_display_qt

        # update vline position
        self.canvas_instance.canvas_annotation_vline.setGeometry(
            self.x, ymax_display_qt, 1, vline_length
        )

        return None
