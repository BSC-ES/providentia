""" Functions and classes to create and format dashboard PyQt elements """

import copy
from difflib import SequenceMatcher
from functools import partial
import platform
import re
from textwrap import wrap

import matplotlib
import numpy as np
from PyQt5 import QtCore, QtWidgets, QtGui
import yaml

from providentia.auxiliar import CURRENT_PATH, join

PROVIDENTIA_ROOT = "/".join(CURRENT_PATH.split("/")[:-1])
# get operating system specific formatting
operating_system = platform.system()
if operating_system == "Darwin":
    formatting_dict = yaml.safe_load(
        open(join(PROVIDENTIA_ROOT, "settings/internal/stylesheet_mac.yaml"))
    )
elif operating_system == "Linux":
    formatting_dict = yaml.safe_load(
        open(join(PROVIDENTIA_ROOT, "settings/internal/stylesheet_linux.yaml"))
    )
elif operating_system in ["Windows", "MINGW32_NT", "MINGW64_NT"]:
    formatting_dict = yaml.safe_load(
        open(join(PROVIDENTIA_ROOT, "settings/internal/stylesheet_windows.yaml"))
    )


def normalise_search_text(text):
    """
    Function which reduces text to the characters a search should compare,
    so that case, spaces, underscores and dashes cannot stop a field being
    found (i.e. "Station Name", "station_name" and "stationname" are equal).

    Parameters
    ----------
    text : str
        Text to normalise

    Returns
    -------
    str
        Text as lowercase letters and digits only
    """

    return re.sub(r"[^0-9a-z]", "", str(text).lower())


def search_match_score(query, label, allow_fuzzy=True):
    """
    Function which scores how well a field label answers a search query,
    returning None when the label is not a match at all.

    Scored rather than simply matched so that results can be ordered by how
    well they answer the query, the closest first. Matching runs through
    progressively looser tests - the whole query as written, then its words in
    any order, then an approximate comparison which tolerates a typo - so an
    exact match is never pushed down the list by a fuzzy one.

    Parameters
    ----------
    query : str
        Text typed into the search box
    label : str
        Field label to test against
    allow_fuzzy : bool, optional
        Whether to fall as far as the approximate, typo-tolerant comparison -
        by far the most expensive test here, run once per label. Left False
        by search_field_labels() on its first pass over a list of labels, so
        it is only ever paid for on a second pass, over the whole list again,
        when that first, cheap pass matched nothing at all - see there.

    Returns
    -------
    float or None
        Score, where lower is a better match, or None if the label does not
        match the query
    """

    query_normalised = normalise_search_text(query)
    label_normalised = normalise_search_text(label)

    # an empty query matches everything, leaving the menu as it was
    if not query_normalised:
        return 0.0
    if not label_normalised:
        return None

    if label_normalised == query_normalised:
        return 0.0
    if label_normalised.startswith(query_normalised):
        return 0.1
    if query_normalised in label_normalised:
        # earlier in the label reads as the better match
        return 0.2 + (
            label_normalised.index(query_normalised) / len(label_normalised)
        )

    # every word of the query somewhere in the label, in any order, so that
    # "class area" finds "area_classification"
    query_words = [
        normalise_search_text(word) for word in re.split(r"[\s_-]+", query.strip())
    ]
    query_words = [word for word in query_words if word]
    if len(query_words) > 1 and all(
        word in label_normalised for word in query_words
    ):
        return 0.5

    if not allow_fuzzy:
        return None

    # approximate match, for a mistyped or half-remembered field name. Held
    # back to longer queries, below which nearly everything looks similar to
    # everything else, and only ever used by search_field_labels() when the
    # query matches nothing as written
    if len(query_normalised) >= 5:
        ratio = SequenceMatcher(None, query_normalised, label_normalised).ratio()
        # also compared against the best window of the label the query could
        # sit in, so a typo still finds a long name the query is only part of
        window = len(query_normalised)
        for start in range(max(1, len(label_normalised) - window + 1)):
            ratio = max(
                ratio,
                SequenceMatcher(
                    None, query_normalised, label_normalised[start : start + window]
                ).ratio(),
            )
        if ratio >= 0.75:
            return 1.0 + (1.0 - ratio)

    return None


def search_field_labels(query, labels):
    """
    Function which picks the field labels answering a search query, ordered
    with the closest match first.

    Approximate matches are only ever offered when nothing matches the query
    as written: a typo should find the field it was meant to be, but should
    not pad a good set of results with loosely similar names.

    Run as two passes rather than one for exactly that reason, and so that
    the (by far the most expensive, a per-label loop of its own) approximate
    comparison is only ever paid for on the rare query that needs it - a
    field of many thousands of labels (e.g. every station name) otherwise
    lagged on every keystroke of an ordinary, plainly-matching query, since a
    single scoring pass ran the approximate comparison for every label that
    the cheap tests above it did not already match.

    Parameters
    ----------
    query : str
        Text typed into the search box
    labels : list
        Field labels to search through

    Returns
    -------
    list
        Indices of the matching labels, best match first
    """

    if not normalise_search_text(query):
        return list(range(len(labels)))

    scored = [
        (score, label_ii)
        for label_ii, label in enumerate(labels)
        for score in [search_match_score(query, label, allow_fuzzy=False)]
        if score is not None
    ]

    # nothing matched as written - worth the far more expensive approximate
    # comparison now, over the whole list again, as the only remaining option
    if not scored:
        scored = [
            (score, label_ii)
            for label_ii, label in enumerate(labels)
            for score in [search_match_score(query, label, allow_fuzzy=True)]
            if score is not None
        ]

    return [label_ii for _, label_ii in sorted(scored)]


def set_formatting(
    PyQt5_obj, format, valid_obj=None, disabled=False, extra_arguments={}
):
    """
    Function that takes a PyQt5 object and applies some defined formatting

    Parameters
    ----------
    PyQt5_obj : object
        PyQt5 element
    format : dict
        Format dictionary
    valid_obj : list
        PyQt5 element to format
    disabled : bool
        Whether we want to format the disabled version of the object
    extra_arguments : list
        Extra arguments

    Returns
    -------
    object
        PyQt5 element with new format
    """

    # initialise style
    full_defined_style = ""

    # iterate through formatting dictionary and apply defined font modifiers/object formatting values
    for obj_type in format:
        if valid_obj:
            if obj_type not in valid_obj:
                continue

        if len(extra_arguments) > 0:
            if obj_type in extra_arguments:
                cut_extra_arguments = extra_arguments[obj_type]
            else:
                cut_extra_arguments = copy.deepcopy(extra_arguments)
        else:
            cut_extra_arguments = {}

        defined_style = ""

        for format_name, format_val in format[obj_type].items():
            if format_name in cut_extra_arguments:
                format_val = cut_extra_arguments[format_name]
                del cut_extra_arguments[format_name]

            if format_name == "height":
                PyQt5_obj.setFixedHeight(int(format_val))
            elif format_name == "width":
                PyQt5_obj.setFixedWidth(int(format_val))
            elif format_name == "min-height":
                PyQt5_obj.setMinimumHeight(int(format_val))
            elif format_name == "min-width":
                PyQt5_obj.setMinimumWidth(int(format_val))
            elif format_name == "max-height":
                PyQt5_obj.setMaximumHeight(int(format_val))
            elif format_name == "max-width":
                PyQt5_obj.setMaximumWidth(int(format_val))
            else:
                defined_style += "{}: {};".format(format_name, format_val)

        # have remaining extra arguments to add?
        if len(cut_extra_arguments) > 0:
            for format_name, format_val in cut_extra_arguments.items():
                if format_name == "height":
                    PyQt5_obj.setFixedHeight(int(format_val))
                elif format_name == "width":
                    PyQt5_obj.setFixedWidth(int(format_val))
                elif format_name == "min-height":
                    PyQt5_obj.setMinimumHeight(int(format_val))
                elif format_name == "min-width":
                    PyQt5_obj.setMinimumWidth(int(format_val))
                elif format_name == "max-height":
                    PyQt5_obj.setMaximumHeight(int(format_val))
                elif format_name == "max-width":
                    PyQt5_obj.setMaximumWidth(int(format_val))
                else:
                    defined_style += "{}: {};".format(format_name, format_val)

        if disabled:
            defined_style = "{}:disabled {{ {} }} ".format(obj_type, defined_style)
        else:
            defined_style = "{} {{ {} }} ".format(obj_type, defined_style)
        full_defined_style += defined_style

    # apply style sheet
    PyQt5_obj.setStyleSheet(full_defined_style)

    return PyQt5_obj


def wrap_tooltip_text(tooltip_text, max_width, format_type):
    """
    Function which takes the text for a tooltip and wraps it by the screen pixel width.
    It does this by estimating the pixel width of the tooltip text (as formatted),
    and then gets the ratio exceedance over the screen pixel width.
    If there is an exceedance (i.e. > 1), the text is then broken into n max_char pieces
    based on the position of the first exceedance in the text
    (i.e. the part of the text which first exceeds the screen pixel width).

    Parameters
    ----------
    tooltip_text : str
        Tooltip text
    max_width : int
        Maximum width accepted
    format_type : str
        Element type

    Returns
    -------
    str
        Updated tooltip text
    """

    tooltip_label = set_formatting(
        QtWidgets.QLabel(text=tooltip_text),
        formatting_dict[format_type],
        valid_obj=["QToolTip"],
    )
    tooltip_width = (
        tooltip_label.fontMetrics().boundingRect(tooltip_label.text()).width()
    )
    if tooltip_width > max_width:
        ratio = tooltip_width / max_width
        max_char = int(np.floor((len(tooltip_text) / ratio) * 1.0))
        tooltip_text = "\n".join(wrap(tooltip_text, max_char))

    return tooltip_text


def center(window):
    """
    Center window
    Reference: https://wiki.qt.io/How_to_Center_a_Window_on_the_Screen

    Parameters
    ----------
    window : MessageBox
        Message box
    """

    window.setGeometry(
        QtWidgets.QStyle.alignedRect(
            QtCore.Qt.LeftToRight,
            QtCore.Qt.AlignCenter,
            window.size(),
            QtWidgets.qApp.desktop().availableGeometry(),
        )
    )


# qualitative colourmaps are a fixed handful of unordered category colours,
# not a continuum - drawn continuously they read as arbitrary bands, and asked
# for more chunks than they have colours they repeat themselves, leaving two
# value ranges the same colour. Nothing errors, so they are excluded here
_QUALITATIVE_COLOURMAPS = {
    "Accent",
    "Dark2",
    "Paired",
    "Pastel1",
    "Pastel2",
    "Set1",
    "Set2",
    "Set3",
    "tab10",
    "tab20",
    "tab20b",
    "tab20c",
}


def get_valid_colourmaps():
    """
    Get the names of every colourmap registered with matplotlib.

    The "_r" (reversed) variants are included: reversing carries meaning here,
    as a statistic whose best value is its maximum is drawn with the reverse
    of the colourmap used for one whose best value is its minimum (see
    settings/colourmaps.yaml), so they have to be selectable.

    Returns
    -------
    list of str
        Sorted colourmap names, valid to pass to matplotlib.colormaps[name]
        or matplotlib.pyplot.get_cmap().
    """

    return sorted(
        name
        for name in matplotlib.colormaps
        if name not in _QUALITATIVE_COLOURMAPS
        and name.removesuffix("_r") not in _QUALITATIVE_COLOURMAPS
    )


def make_colourmap_icon(name, width=64, height=13):
    """
    Render a small horizontal gradient swatch showing a matplotlib
    colourmap's actual colours, for use as a QComboBox item icon.

    Parameters
    ----------
    name : str
        A valid matplotlib colourmap name (see get_valid_colourmaps()).
    width : int, optional
        Icon width in pixels (default is 64).
    height : int, optional
        Icon height in pixels (default is 13).

    Returns
    -------
    QtGui.QIcon
        Icon showing the colourmap as a rounded gradient swatch.
    """

    cmap = matplotlib.colormaps[name]
    pixmap = QtGui.QPixmap(width, height)
    pixmap.fill(QtCore.Qt.transparent)

    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.Antialiasing)
    gradient = QtGui.QLinearGradient(0, 0, width, 0)
    n_stops = 16
    for stop_ii in range(n_stops):
        fraction = stop_ii / (n_stops - 1)
        r, g, b, a = cmap(fraction)
        gradient.setColorAt(fraction, QtGui.QColor.fromRgbF(r, g, b, a))
    painter.setPen(QtGui.QPen(QtGui.QColor("#C0CBD1"), 1))
    painter.setBrush(QtGui.QBrush(gradient))
    painter.drawRoundedRect(0, 0, width - 1, height - 1, 3, 3)
    painter.end()

    return QtGui.QIcon(pixmap)


def select_colourmap(combobox, name):
    """
    Select a colourmap in a combobox already populated by
    populate_colourmap_combobox(), leaving its items alone.

    Rebuilding the list to change the selection means regenerating every
    swatch icon, which is wasted work when only the selection is changing.

    Parameters
    ----------
    combobox : QtWidgets.QComboBox
        Combobox to change the selection of.
    name : str
        Colourmap name to select. Ignored if it is not in the combobox.
    """

    index = combobox.findText(name)
    if index != -1:
        combobox.setCurrentIndex(index)


def populate_colourmap_combobox(combobox, current=None):
    """
    Fill a QComboBox with every valid matplotlib colourmap, each shown with a
    gradient swatch icon of its actual colours.

    Parameters
    ----------
    combobox : QtWidgets.QComboBox
        Combobox to populate (cleared first).
    current : str, optional
        Colourmap name to select initially. Falls back to the first (alpha-
        betically) colourmap if this isn't a valid colourmap name.
    """

    combobox.clear()
    combobox.setIconSize(QtCore.QSize(64, 13))
    valid_colourmaps = get_valid_colourmaps()
    for name in valid_colourmaps:
        combobox.addItem(make_colourmap_icon(name), name)

    # ComboBox is editable, and setCurrentText() sets its internal line edit
    # without syncing the closed-box icon, leaving a stale swatch -
    # setCurrentIndex() goes through the real selection path instead
    target_text = current if current in valid_colourmaps else valid_colourmaps[0]
    combobox.setCurrentIndex(valid_colourmaps.index(target_text))


# cartopy projections that construct with no required arguments, render a
# sensible whole-world view by default, and aren't a near-duplicate or a
# regional CRS. Checked by rendering each candidate, not by name alone
VALID_PROJECTIONS = [
    "Aitoff",
    "AzimuthalEquidistant",
    "EckertI",
    "EckertII",
    "EckertIII",
    "EckertIV",
    "EckertV",
    "EckertVI",
    "EqualEarth",
    "Hammer",
    "InterruptedGoodeHomolosine",
    "LambertAzimuthalEqualArea",
    "LambertCylindrical",
    "Mercator",
    "Miller",
    "Mollweide",
    "NorthPolarStereo",
    "Orthographic",
    "PlateCarree",
    "Robinson",
    "Sinusoidal",
    "SouthPolarStereo",
    "Stereographic",
]


def get_valid_projections():
    """
    Get the names of the cartopy (cartopy.crs) map projections offered in the
    map settings menu - see the VALID_PROJECTIONS comment for how this list
    was curated.

    Returns
    -------
    list of str
        Sorted projection names, valid to pass to getattr(cartopy.crs, name).
    """

    return sorted(VALID_PROJECTIONS)


def make_projection_icon(name, width=34, height=18):
    """
    Build a small icon showing a projection's actual world outline (land +
    the projection's own boundary shape - an ellipse for Mollweide, a circle
    for Orthographic, a hexagon for EckertI, etc - pre-rendered with cartopy;
    see providentia/resources/projections/), for use as a QComboBox item
    icon. Unlike make_colourmap_icon()/make_colour_icon(), there's no
    separate card frame drawn here - the projection's boundary *is* the
    icon's shape, which is the whole point of showing it.

    Parameters
    ----------
    name : str
        A valid projection name (see get_valid_projections()).
    width : int, optional
        Icon width in pixels (default is 34).
    height : int, optional
        Icon height in pixels (default is 18).

    Returns
    -------
    QtGui.QIcon
        Icon showing the projection's true outline shape.
    """

    pixmap = QtGui.QPixmap(width, height)
    pixmap.fill(QtCore.Qt.transparent)

    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.Antialiasing)
    painter.setRenderHint(QtGui.QPainter.SmoothPixmapTransform)

    projection_path = join(
        PROVIDENTIA_ROOT, f"providentia/resources/projections/{name}.png"
    )
    projection_pixmap = QtGui.QPixmap(projection_path)
    if not projection_pixmap.isNull():
        projection_pixmap = projection_pixmap.scaled(
            width,
            height,
            QtCore.Qt.KeepAspectRatio,
            QtCore.Qt.SmoothTransformation,
        )
        x = (width - projection_pixmap.width()) // 2
        y = (height - projection_pixmap.height()) // 2
        painter.drawPixmap(x, y, projection_pixmap)
    painter.end()

    return QtGui.QIcon(pixmap)


def populate_projection_combobox(combobox, current):
    """
    Fill a QComboBox with every valid map projection, each shown with an icon
    of its actual world outline.

    Parameters
    ----------
    combobox : QtWidgets.QComboBox
        Combobox to populate (cleared first).
    current : str
        Projection name to select initially. Falls back to the first
        (alphabetically) projection if this isn't a valid projection name.
    """

    combobox.clear()
    combobox.setIconSize(QtCore.QSize(34, 18))
    valid_projections = get_valid_projections()
    for name in valid_projections:
        combobox.addItem(make_projection_icon(name), name)

    # see populate_colourmap_combobox() for why this is setCurrentIndex(),
    # not setCurrentText()
    target_text = current if current in valid_projections else valid_projections[0]
    combobox.setCurrentIndex(valid_projections.index(target_text))


# curated so every option reads clearly at swatch size and stays muted enough
# not to compete with the station colours on top. Hex matches the
# map.land_polygon/map.ocean_polygon defaults in plot_characteristics.yaml, so
# the shipped colours show as selected on first open
LAND_COLOUR_OPTIONS = {
    "Light grey": "#D9D9D9",
    "Grey": "#B0B0B0",
    "Beige": "#E8DCC8",
    "Sand": "#EDE0C8",
    "Light green": "#C8DCC0",
    "Green": "#8FBC8F",
    "Off white": "#E8E8E8",
    "White": "#FFFFFF",
    "Charcoal": "#5A5A5A",
    # dark end of the range, so the map as a whole can be turned dark
    # (pair with one of the dark ocean options below) rather than only
    # ever being a light basemap with a single dark shade available
    "Dark grey": "#3A3A3A",
    "Dark slate": "#2F3640",
    "Dark green": "#2B3B30",
    "Near black": "#1E1E1E",
}

OCEAN_COLOUR_OPTIONS = {
    "Soft blue": "#DCE6ED",
    "Sky blue": "#AED6F1",
    "Steel blue": "#7FA8C9",
    "Deep blue": "#2C5F8A",
    "Navy": "#1B4B8F",
    "Teal": "#4FB6CE",
    "Grey": "#D5D8DC",
    "Off white": "#FAFAFA",
    "White": "#FFFFFF",
    # dark counterparts to the land options above
    "Dark grey": "#33383D",
    "Dark slate": "#22303C",
    "Midnight": "#121C26",
    "Near black": "#0D0D0D",
}


# ready-made land/ocean/colourmap combinations, chosen so the basemap stays
# quiet enough for the station colours to carry the information. "Custom" is
# not a preset - it is what the selector shows once land, ocean or colourmap


def make_colour_icon(hex_colour, width=34, height=13):
    """
    Render a small flat colour swatch icon, for use as a QComboBox item icon
    (see make_colourmap_icon()/make_projection_icon() for the same idea
    applied to colourmaps and projections).

    Parameters
    ----------
    hex_colour : str
        Colour as a "#RRGGBB" string.
    width : int, optional
        Icon width in pixels (default is 34).
    height : int, optional
        Icon height in pixels (default is 13).

    Returns
    -------
    QtGui.QIcon
        Icon showing the colour as a rounded swatch.
    """

    pixmap = QtGui.QPixmap(width, height)
    pixmap.fill(QtCore.Qt.transparent)

    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.Antialiasing)
    painter.setPen(QtGui.QPen(QtGui.QColor("#C0CBD1"), 1))
    painter.setBrush(QtGui.QBrush(QtGui.QColor(hex_colour)))
    painter.drawRoundedRect(0, 0, width - 1, height - 1, 3, 3)
    painter.end()

    return QtGui.QIcon(pixmap)


def populate_colour_combobox(combobox, options, current):
    """
    Fill a QComboBox with a curated set of named colours, each shown with a
    flat swatch icon.

    Parameters
    ----------
    combobox : QtWidgets.QComboBox
        Combobox to populate (cleared first).
    options : dict
        {label: "#RRGGBB"} - e.g. LAND_COLOUR_OPTIONS or OCEAN_COLOUR_OPTIONS.
    current : str
        Hex colour to select initially (matched against options' values, not
        the labels). Falls back to the first option if not a match.
    """

    combobox.clear()
    combobox.setIconSize(QtCore.QSize(34, 13))
    labels = list(options.keys())
    for label, hex_colour in options.items():
        combobox.addItem(make_colour_icon(hex_colour), label)

    matches = [label for label, hex_colour in options.items() if hex_colour == current]
    target_label = matches[0] if matches else labels[0]
    # see populate_colourmap_combobox() for why this is setCurrentIndex(),
    # not setCurrentText()
    combobox.setCurrentIndex(labels.index(target_label))


class MenuLineEdit(QtWidgets.QLineEdit):
    """
    Settings-menu line edit that reliably reports when the user has
    finished with it - on Enter, and on clicking away.

    Qt's own editingFinished() is not dependable for that second case
    here: it is suppressed whenever a validator judges the current text
    merely "intermediate" rather than acceptable, which an empty integer
    field always is. Clearing one of the colourbar fields back to
    automatic and clicking elsewhere therefore emitted nothing at all,
    and the map kept the old value until Enter was pressed in the field.

    Emits at most once per actual change: the text at the last commit is
    remembered, so tabbing or clicking through a field without touching
    it doesn't trigger a redraw, and pressing Enter and then clicking
    away doesn't trigger two.
    """

    committed = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        """
        Initialise class

        Parameters
        ----------
        parent : QtWidgets.QWidget, optional
            Parent widget.
        """

        super().__init__(parent)
        self._committed_text = self.text()
        self.returnPressed.connect(self._commit)

    def _commit(self):
        """Emit committed() if the text has changed since the last one."""

        if self.text() != self._committed_text:
            self._committed_text = self.text()
            self.committed.emit()

    def setText(self, text):
        """
        Set the text, treating it as the new committed baseline.

        Values written programmatically (e.g. the resolved colourbar
        limits written back after every redraw) are not user edits, so
        they must not count as a pending change that a later focus-out
        would then re-emit.

        Parameters
        ----------
        text : str
            Text to set.
        """

        super().setText(text)
        self._committed_text = self.text()

    def clear(self):
        """Clear the text, treating the empty value as committed."""

        super().clear()
        self._committed_text = self.text()

    def has_pending_edit(self):
        """
        Whether the text has been changed since it was last committed.

        Returns
        -------
        bool
            True if there is an uncommitted edit.
        """

        return self.text() != self._committed_text

    def mark_committed(self):
        """
        Record the current text as committed, without emitting anything.

        For callers that apply the value themselves rather than through
        this widget's signal - see Canvas.commit_map_pending_edits().
        """

        self._committed_text = self.text()

        return None

    def commit_pending(self):
        """
        Commit straight away if the text has changed since the last
        commit, rather than waiting for Enter or focus to move.

        For the cases neither of those covers - the panel holding the
        field being closed while it still has focus - see
        Canvas.commit_map_pending_edits(). Runs inline rather than
        deferred: there is no click still in flight to get out of the way
        of, and the caller is about to hide the widget.
        """

        self._commit()

        return None

    def focusOutEvent(self, event):
        """
        Commit on losing focus, then hand on to the default handler.

        The commit is deferred to the next pass of the event loop rather
        than run inline. Focus is lost *during* the mouse press of
        whatever was clicked next, and committing here redraws the map
        synchronously - which swallowed the rest of that click, so a
        button next to one of these fields (the colourbar resets) had to
        be pressed twice: once to commit the field, again to actually
        activate. Deferring lets the click complete first.

        An application-wide event filter that committed on the press
        itself was tried instead, and was worse for the same reason: the
        handler it invoked re-entered the event loop (processEvents) in
        the middle of Qt delivering that press, so the button never saw a
        complete press/release pair and its clicked signal never fired at
        all. Nothing may do synchronous work during delivery - the value
        is applied afterwards here, and, for menu controls, ahead of the
        handler by SettingsMenu.connect().

        Parameters
        ----------
        event : QtGui.QFocusEvent
            Focus event.
        """

        if self.text() != self._committed_text:
            QtCore.QTimer.singleShot(0, self._commit)
        super().focusOutEvent(event)


class ComboBox(QtWidgets.QComboBox):
    """Modify default class of PyQT5 combobox."""

    def __init__(self, parent=None):
        """
        Initialise class

        Parameters
        ----------
        parent : object
            Dashboard instance
        """

        super(ComboBox, self).__init__(parent)

        # setMaxVisibleItems only works if the box is editable
        # this creates a line edit that we need to overwrite
        self.setEditable(True)
        # self.AdjustToContents
        self.setSizeAdjustPolicy(self.AdjustToMinimumContentsLengthWithIcon)
        self.setMaxVisibleItems(20)
        self.AdjustToContents

        # overwrite default line edit by an invisible one
        self.lineEdit().setFrame(False)
        self.lineEdit().setReadOnly(True)
        self.currentTextChanged.connect(self.fixCursorPosition)

    def fixCursorPosition(self):
        """
        Move (invisible) cursor to first position, so that a name too long for
        the box is shown from its start rather than its end.

        A line edit shows the text around wherever its cursor sits, which
        leaves a long name identifiable only by its ending. Applied whenever
        the text or the width changes, as either can leave the box scrolled.
        """

        self.lineEdit().setCursorPosition(0)

    def setCurrentText(self, text):
        """
        Set the text shown by the box

        Parameters
        ----------
        text : str
            Text to show
        """

        # setting the text leaves the cursor at its end, and says nothing when
        # the text has not actually changed, so the box is put back to its
        # start here rather than left to currentTextChanged
        super().setCurrentText(text)
        self.fixCursorPosition()

    def setCurrentIndex(self, index):
        """
        Set which of the box's options is shown

        Parameters
        ----------
        index : int
            Index of the option
        """

        super().setCurrentIndex(index)
        self.fixCursorPosition()

    def resizeEvent(self, event):
        """
        Handle the box being resized

        Parameters
        ----------
        event : QResizeEvent
            Resize event
        """

        super().resizeEvent(event)
        self.fixCursorPosition()

    def showEvent(self, event):
        """
        Handle the box being shown

        Parameters
        ----------
        event : QShowEvent
            Show event
        """

        super().showEvent(event)
        self.fixCursorPosition()

    def showPopup(self):
        """
        Show pop-up
        """

        # set index of selected choice to highlight it
        text = self.lineEdit().text()
        index = self.findText(text, QtCore.Qt.MatchFixedString)
        self.setCurrentIndex(index)

        # show pop-up
        super().showPopup()

        # increase the width of the elements on popup so they can be read
        self.view().setMinimumWidth(self.view().sizeHintForColumn(0) + 10)

        # add vertical scroll bar
        self.view().setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)


class CheckableComboBox(QtWidgets.QComboBox):
    def __init__(self, *args, **kwargs):
        """
        Initialise class
        """

        super().__init__(*args, **kwargs)

        # make the combo editable to set a custom text, but readonly
        self.setEditable(True)
        self.lineEdit().setReadOnly(True)
        self.lineEdit().setPlaceholderText("Select option/s:")
        self.setMaxVisibleItems(20)
        self.currentTextChanged.connect(self.fixCursorPosition)

        # make the lineedit the same color as QComboBox
        palette = QtWidgets.QApplication.palette()
        palette.setBrush(QtGui.QPalette.Base, palette.button())
        self.lineEdit().setPalette(palette)

        # update the text when an item is toggled
        self.model().dataChanged.connect(self.updateText)

        # hide and show popup when clicking the line edit
        self.lineEdit().installEventFilter(self)
        self.closeOnLineEditClick = False

        # prevent popup from closing when clicking on an item
        self.view().viewport().installEventFilter(self)

    def fixCursorPosition(self):
        """
        Move (invisible) cursor to first position, so that a name too long for
        the box is shown from its start rather than its end.

        A line edit shows the text around wherever its cursor sits, which
        leaves a long name identifiable only by its ending. Applied whenever
        the text or the width changes, as either can leave the box scrolled.
        """

        self.lineEdit().setCursorPosition(0)

    def resizeEvent(self, event):
        """
        Resize event after updating text

        Parameters
        ----------
        event : QResizeEvent
            Resize event
        """

        # recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)
        self.fixCursorPosition()

    def eventFilter(self, obj, event):
        """
        Custom multi-select dropdown where clicking items checks/unchecks them while keeping the
        popup visually smooth.

        Parameters
        ----------
        obj : object
            Element
        event : QtCore.QEvent
            Event
        """

        if obj == self.lineEdit():
            if event.type() == QtCore.QEvent.MouseButtonRelease:
                if self.closeOnLineEditClick:
                    self.hidePopup()
                else:
                    self.showPopup()
                return True
            return False

        if obj == self.view().viewport():
            if event.type() == QtCore.QEvent.MouseButtonRelease:
                index = self.view().indexAt(event.pos())
                if not index.isValid():
                    return False

                item = self.model().item(index.row())
                if not (item.flags() & QtCore.Qt.ItemIsEnabled):
                    return True

                # Toggle check state
                item.setCheckState(
                    QtCore.Qt.Unchecked
                    if item.checkState() == QtCore.Qt.Checked
                    else QtCore.Qt.Checked
                )

                # Save current visual state
                row_to_restore = index.row()
                scroll_value = self.view().verticalScrollBar().value()

                # Hide now; reopen after the layout/menu has settled
                self.hidePopup()

                def _reopen_and_restore():
                    # Block repaints during reopen to avoid flicker
                    self.setUpdatesEnabled(False)
                    self.view().setUpdatesEnabled(False)

                    self.showPopup()

                    # Restore scroll & selection
                    try:
                        self.view().verticalScrollBar().setValue(scroll_value)
                    except RuntimeError:
                        pass  # view may be recreated

                    if 0 <= row_to_restore < self.model().rowCount():
                        idx = self.model().index(row_to_restore, 0)
                        if idx.isValid():
                            self.view().setCurrentIndex(idx)
                            self.view().scrollTo(idx)

                    self.view().setUpdatesEnabled(True)
                    self.setUpdatesEnabled(True)

                # Defer to next cycle so geometry changes (from text update) are applied
                QtCore.QTimer.singleShot(0, _reopen_and_restore)
                return True

        return False

    def showPopup(self):
        """
        Custom show pop up
        """

        model = self.model()

        # Find the first enabled item
        first_enabled_index = -1
        for i in range(model.rowCount()):
            item = model.item(i)
            if item.flags() & QtCore.Qt.ItemIsEnabled:
                first_enabled_index = i
                break

        # Open the popup first
        super().showPopup()

        # Highlight the first enabled item in the popup view only
        if first_enabled_index != -1:
            idx = model.index(first_enabled_index, 0)
            self.view().selectionModel().clearSelection()  # clear any existing selection
            self.view().selectionModel().setCurrentIndex(
                idx, QtCore.QItemSelectionModel.SelectCurrent
            )
            self.view().scrollTo(idx)

        # Adjust popup width
        width = self.view().sizeHintForColumn(0) + 40
        self.view().setMinimumWidth(width)
        self.closeOnLineEditClick = True

    def hidePopup(self):
        """
        Custom hide pop up
        """

        super().hidePopup()
        self.startTimer(100)
        self.updateText()

    def timerEvent(self, event):
        """
        Stop timer and disable closing the popup

        Parameters
        ----------
        event : QtCore.QEvent
            Event
        """

        self.killTimer(event.timerId())
        self.closeOnLineEditClick = False

    def updateText(self):
        """
        Show checked elements as one string in line
        """

        texts = [
            self.model().item(i).text()
            for i in range(self.model().rowCount())
            if self.model().item(i).checkState() == QtCore.Qt.Checked
        ]
        self.lineEdit().setText(", ".join(texts))

    def addItem(self, text, data=None, enabled=True):
        """
        Add a single checkable item to the model

        Parameters
        ----------
        text : str
            Display text of the item
        data : int, None
            Associated data
        enabled : bool, default True
            If False, item is disabled and grayed out.
        """

        item = QtGui.QStandardItem()
        item.setText(text)
        item.setData(data if data is not None else text)
        flags = QtCore.Qt.ItemIsUserCheckable
        if enabled:
            flags |= QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
        item.setFlags(flags)
        item.setData(QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)

        if not enabled:
            item.setData(
                QtGui.QBrush(QtGui.QColor(150, 150, 150)), QtCore.Qt.ForegroundRole
            )

        self.model().appendRow(item)

        # Highlight the first enabled item if nothing is selected
        if self.currentIndex() == -1:
            for i in range(self.model().rowCount()):
                first_item = self.model().item(i)
                if first_item.flags() & QtCore.Qt.ItemIsEnabled:
                    self.setCurrentIndex(i)
                    break

    def addItems(self, texts, datalist=None, enabled_list=None):
        """
        Add multiple checkable items to the model

        Parameters
        ----------
        texts : list of str
            Display texts for the items
        datalist : list
            Data for each item
        enabled_list : list
            Enabled state for each item
        """

        for i, text in enumerate(texts):
            data = datalist[i] if datalist and i < len(datalist) else None
            enabled = (
                enabled_list[i] if enabled_list and i < len(enabled_list) else True
            )
            self.addItem(text, data, enabled)

    def currentData(self, all=False):
        """
        Return item data from the model

        Parameters
        ----------
        all : bool
            If True, return data for all items
            If False, return data only for checked items

        Returns
        -------
        list
            List of item data values
        """

        return [
            self.model().item(i).data()
            for i in range(self.model().rowCount())
            if all or self.model().item(i).checkState() == QtCore.Qt.Checked
        ]


class QVLine(QtWidgets.QFrame):
    """
    Define class that generates vertical separator line
    """

    def __init__(self, parent=None):
        super(QVLine, self).__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class LegendInlineEditor(QtWidgets.QLineEdit):
    """
    Small QLineEdit overlaid directly on top of a matplotlib legend label
    for in-place renaming (see rename_legend_label() in
    dashboard_interactivity.py), instead of a separate pop-up dialog.

    QLineEdit already emits editingFinished on both Enter and losing focus,
    which covers "commit"; it has no equivalent for "cancel", so this adds
    an escapePressed signal for that.
    """

    escapePressed = QtCore.pyqtSignal()

    def keyPressEvent(self, event):
        """
        Emit escapePressed on the Escape key, otherwise behave as normal.

        Parameters
        ----------
        event : QKeyEvent
            The key press event.
        """

        if event.key() == QtCore.Qt.Key_Escape:
            self.escapePressed.emit()
        else:
            super(LegendInlineEditor, self).keyPressEvent(event)


class MenuEditCommitFilter(QtCore.QObject):
    """
    Installed once, application-wide, so that clicking anywhere outside
    the settings field being edited applies its value.

    Focus alone is not a dependable signal for this. Whether a click moves
    focus at all depends on the widget it lands on and its focus policy -
    plain panel background, a label, the canvas and a button all behave
    differently, and on macOS several of them do not take focus at all -
    so "clicking on idle space" reached the field's focus-out handler in
    some places and not others.

    This watches the mouse *release* rather than the press, and only
    schedules the commit rather than running it. Both matter: an earlier
    version committed during the press, and the redraw that triggered ran
    in the middle of Qt delivering the click, so the button underneath
    never completed its press/release pair. By release the click is
    already done, and deferring to the next pass of the event loop keeps
    any redraw out of event delivery entirely.
    """

    def __init__(self, parent=None):
        """
        Initialise class

        Parameters
        ----------
        parent : QtCore.QObject, optional
            Parent object.
        """

        super(MenuEditCommitFilter, self).__init__(parent)

    def eventFilter(self, obj, event):
        """
        Schedule a commit of the focused settings field when a click
        finishes anywhere else.

        Parameters
        ----------
        obj : QtCore.QObject
            Object the event was sent to.
        event : QtCore.QEvent
            The event.

        Returns
        -------
        bool
            Always False - this only observes, it never consumes.
        """

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            focused = QtWidgets.QApplication.focusWidget()
            if (
                isinstance(focused, MenuLineEdit)
                and focused is not obj
                and focused.has_pending_edit()
                and not (
                    isinstance(obj, QtWidgets.QWidget) and focused.isAncestorOf(obj)
                )
            ):
                QtCore.QTimer.singleShot(0, focused.commit_pending)

        return False


class LegendEditorCommitFilter(QtCore.QObject):
    """
    Installed application-wide for the lifetime of a LegendInlineEditor -
    a click landing anywhere other than the editor itself commits the
    rename, the same as pressing Enter. Needed because the editor sits
    over a matplotlib canvas: clicking elsewhere on that canvas is handled
    by matplotlib's own pick/event machinery rather than normal Qt
    click-to-focus, so the editor doesn't reliably lose focus (and
    therefore never emits editingFinished) just from clicking away.
    """

    def __init__(self, editor, on_outside_press):
        """Initialise class

        Parameters
        ----------
        editor : LegendInlineEditor
            The editor this filter is guarding - a press delivered to any
            other widget counts as "outside".
        on_outside_press : callable
            Called (no arguments) the first time a press outside the
            editor is observed.
        """

        super(LegendEditorCommitFilter, self).__init__(editor)
        self.editor = editor
        self.on_outside_press = on_outside_press

    def eventFilter(self, obj, event):
        """
        Watch every mouse press/double-click application-wide, triggering the
        commit callback for anything not delivered to the editor itself.

        Parameters
        ----------
        obj : QtCore.QObject
            Object the event was delivered to
        event : QtCore.QEvent
            Event being filtered

        Returns
        -------
        bool
            Always False, as this only observes and never consumes the event
        """

        if event.type() in (QtCore.QEvent.MouseButtonPress, QtCore.QEvent.MouseButtonDblClick):
            if obj is not self.editor:
                self.on_outside_press()
        return False


class Switch(QtWidgets.QPushButton):
    """Define class that generates switch buttons."""

    def __init__(self, parent=None):
        """
        Initialise class

        Parameters
        ----------
        parent : object
            Dashboard instance
        """

        super(Switch, self).__init__(parent)
        self.setCheckable(True)

    def paintEvent(self, event):
        """
        Define switch properties

        Parameters
        ----------
        event : QtCore.QEvent
            Event
        """

        # set switch properties
        radius = 9
        width = 20
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.translate(self.rect().center())

        # add grey border to switch main box
        painter.setPen(QtGui.QPen(QtCore.Qt.gray))

        # set white background
        painter.setBrush(QtCore.Qt.white)
        painter.drawRoundedRect(
            QtCore.QRect(-width, -radius, 2 * width, 2 * radius), radius, radius
        )

        # set colours and labels on switch
        label = "ON" if self.isChecked() else "OFF"
        bg_colour = QtCore.Qt.black if self.isChecked() else QtCore.Qt.gray
        text_colour = QtCore.Qt.white if self.isChecked() else QtCore.Qt.black

        # set switch background color
        painter.setBrush(QtGui.QBrush(bg_colour))

        # remove switch border color
        painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))

        # change position depending on check
        sw_rect = QtCore.QRect(-radius, -radius, width + radius, 2 * radius)
        if not self.isChecked():
            sw_rect.moveLeft(-width)
        painter.drawRoundedRect(sw_rect, radius, radius)

        # add label (ON / OFF)
        painter.setPen(QtGui.QPen(text_colour))
        painter.drawText(sw_rect, QtCore.Qt.AlignCenter, label)


class MessageBox(QtWidgets.QWidget):
    def __init__(self, msg, parent=None, confirmation=False):
        """
        Initialise class

        Parameters
        ----------
        msg : str
            Text on message box
        parent : object
            Dashboard instance
        confirmation : bool
            Indicates whether we want to ask a question with Yes or No as an answer
        """

        super().__init__(parent)

        self.result = False
        msg_box = self.create_msg_box(msg, confirmation)
        if msg_box is not None:
            layout = QtWidgets.QVBoxLayout(self)
            layout.addWidget(msg_box)
            center(self)

    def create_msg_box(self, msg, confirmation=False):
        """
        Create message box

        Parameters
        ----------
        msg : str
            Text on message box
        confirmation : bool
            Indicates whether we want to ask a question with Yes or No as an answer
        """

        # add warning box
        msg_box = QtWidgets.QMessageBox()
        msg_box.setWindowTitle("Warning")
        msg_box.setText(msg)

        if confirmation:
            # add yes button
            yes_button = set_formatting(
                QtWidgets.QPushButton("Yes"), formatting_dict["popup_button"]
            )

            # add no button
            no_button = set_formatting(
                QtWidgets.QPushButton("No"), formatting_dict["popup_button"]
            )

            msg_box.addButton(yes_button, QtWidgets.QMessageBox.AcceptRole)
            msg_box.addButton(no_button, QtWidgets.QMessageBox.RejectRole)

        else:
            # add ok button
            ok_button = set_formatting(
                QtWidgets.QPushButton("OK"), formatting_dict["popup_button"]
            )

            msg_box.addButton(ok_button, QtWidgets.QMessageBox.AcceptRole)

        # create wrapper to center
        wrapper = partial(center, msg_box)
        QtCore.QTimer.singleShot(0, wrapper)

        result = msg_box.exec_()

        if confirmation:
            self.result = result == QtWidgets.QMessageBox.AcceptRole


class InputDialog(QtWidgets.QWidget):
    def __init__(self, read_instance, title, msg, options, parent=None):
        """
        Initialise class

        Parameters
        ----------
        read_instance : object
            Instance of class Dashboard or Report
        title : str
            Title of input dialog
        msg : str
            Text of input dialog
        options : list
            Dialog options
        parent : object
            Dashboard instance
        """

        super().__init__(parent)

        dialog = self.create_dialog_box(read_instance, title, msg, options)
        if dialog is not None:
            layout = QtWidgets.QVBoxLayout(self)
            layout.addWidget(dialog)

    def create_dialog_box(self, read_instance, title, msg, options):
        dialog = QtWidgets.QInputDialog(self)
        self.selected_option, self.okpressed = dialog.getItem(
            read_instance, title, msg, options, 0, False
        )
        if not self.okpressed:
            return


class CheckDialog(QtWidgets.QDialog):
    def __init__(self, items):
        super().__init__()

        self.setWindowTitle("Select Items")

        layout = QtWidgets.QVBoxLayout(self)

        self.list_widget = QtWidgets.QListWidget()

        # remove focus
        self.list_widget.setFocusPolicy(QtCore.Qt.NoFocus)

        for item_text in items:
            item = QtWidgets.QListWidgetItem(item_text)

            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.Unchecked)

            self.list_widget.addItem(item)

        self.list_widget.itemClicked.connect(self.toggle_item)

        layout.addWidget(self.list_widget)

        button = QtWidgets.QPushButton("OK")
        button.clicked.connect(self.accept)
        layout.addWidget(button)

    def toggle_item(self, item):
        item.setCheckState(
            QtCore.Qt.Unchecked
            if item.checkState() == QtCore.Qt.Checked
            else QtCore.Qt.Checked
        )

    def get_checked_items(self):
        checked = []

        for i in range(self.list_widget.count()):
            item = self.list_widget.item(i)

            if item.checkState() == QtCore.Qt.Checked:
                checked.append(item.text())

        return checked


def create_custom_cursor(size=24):
    """
    Create custom loading cursor (Providentia logo).

    Parameters
    ----------
    size : int
        Size of the cursor in pixels
    """

    pix = QtGui.QPixmap(join(PROVIDENTIA_ROOT, "assets/logo.png"))
    pix = pix.scaled(
        size, size, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation
    )
    return QtGui.QCursor(pix)


def set_cursor(cursor_function, function, size=24):
    """
    Set custom cursor when performing a function that takes time to execute,
    and restore it to normal when finished.

    Parameters
    ----------
    cursor_function : str
        Name of the function which currently owns the cursor, to avoid resetting the cursor if it is already set by another function
    function : str
        Name of the function which is trying to set the cursor
    size : int
        Size of the cursor in pixels

    Returns
    -------
    str
        Function which owns thr cursor, to be used for restoring the cursor to normal when function is finished
    """

    cursor = QtWidgets.QApplication.overrideCursor()
    if cursor is None:
        pix = create_custom_cursor(size)
        QtWidgets.QApplication.setOverrideCursor(pix)
        return function
    else:
        return cursor_function


def unset_cursor(cursor_function, function):
    """
    Unset custom cursor if the function which owns the cursor is the one trying to restore it,
    to avoid restoring the cursor to normal if it is still being used by another function.

    Parameters
    ----------
    cursor_function : str
        Name of the function which currently owns the cursor, to avoid resetting the cursor if it is already set by another function
    function : str
        Name of the function which is trying to set the cursor
    """
    if cursor_function == function:
        QtWidgets.QApplication.restoreOverrideCursor()
