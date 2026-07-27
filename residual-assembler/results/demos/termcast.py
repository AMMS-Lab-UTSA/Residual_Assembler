#!/usr/bin/env python3
"""termcast -- a scripted terminal-screencast renderer (PIL frames -> ffmpeg mp4).

Renders a macOS-style terminal window (title bar + traffic lights) and animates a
scripted session: green ``# ...`` narration, bold commands typed character by
character, and colorized output. Produces an .mp4 with no external service --
frames are drawn with PIL and piped straight into the ffmpeg binary bundled with
imageio_ffmpeg (or $IMAGEIO_FFMPEG_EXE / a system ffmpeg).

Used by demo_program1.py (umat_oti) and demo_program2.py (resasm) to make the
two demo videos in one consistent style.

Span model: a "line" is a list of (text, color_key[, bold]) tuples. color_key
indexes PALETTE. Bold defaults False; the 'cmd' key renders bold automatically.
"""
import os
import shutil
import subprocess

from PIL import Image, ImageDraw, ImageFont

# --------------------------------------------------------------------------- palette (sampled from the reference screencast)
PALETTE = {
    "bg":      (24, 23, 28),
    "bar":     (40, 38, 46),
    "fg":      (208, 210, 205),   # default text
    "gray":    (110, 112, 122),   # line numbers / title / dim
    "comment": (129, 178, 120),   # green "# ..." narration
    "green":   (124, 197, 125),   # success / PASS
    "cmd":     (233, 236, 231),   # bold command (bright white)
    "puser":   (95, 200, 175),    # user@host (teal)
    "ppath":   (101, 160, 247),   # :~/path (blue)
    "psym":    (150, 155, 165),   # $
    "blue":    (136, 177, 243),   # file names
    "cyan":    (99, 194, 204),    # role var (STRESS / stress)
    "orange":  (223, 148, 78),    # role var (DSTRAN / seed)
    "purple":  (191, 141, 201),   # role var (DDSDDE / tangent)
    "yellow":  (226, 196, 98),    # numbers / highlight
    "red":     (226, 112, 102),
    "key":     (130, 170, 235),   # json keys
}
_LIGHTS = [(255, 93, 85), (254, 188, 44), (38, 201, 64)]
_BOLD_KEYS = {"cmd"}


def _ffmpeg():
    for c in (os.environ.get("IMAGEIO_FFMPEG_EXE"),
              "/home/ammslab3/.local/lib/python3.9/site-packages/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v7.0.2",
              shutil.which("ffmpeg")):
        if c and os.path.exists(c):
            return c
    raise RuntimeError("no ffmpeg binary found (set IMAGEIO_FFMPEG_EXE)")


def _mono():
    for c in ("/opt/intel/oneapi/intelpython/latest/lib/python3.9/site-packages/matplotlib/mpl-data/fonts/ttf/%s",
              "/usr/share/fonts/truetype/dejavu/%s"):
        base = c % "DejaVuSansMono.ttf"
        if os.path.exists(base):
            return base, (c % "DejaVuSansMono-Bold.ttf")
    import matplotlib
    d = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data/fonts/ttf")
    return os.path.join(d, "DejaVuSansMono.ttf"), os.path.join(d, "DejaVuSansMono-Bold.ttf")


class Termcast:
    def __init__(self, out_path, title, w=1180, h=720, fontsize=20, fps=24, crf=20):
        self.W, self.H, self.fps = w, h, fps
        reg, bold = _mono()
        self.font = ImageFont.truetype(reg, fontsize)
        self.fontb = ImageFont.truetype(bold, fontsize)
        self.tfont = ImageFont.truetype(reg, 17)
        self.cw = self.font.getlength("M")
        self.lh = int(round(fontsize * 1.33))
        self.pad_x = 26
        self.bar_h = 44
        self.top = self.bar_h + 16
        self.title = title
        self.lines = []                       # committed lines (each: list of spans)
        self.max_lines = (h - self.top - 8) // self.lh
        self.proc = subprocess.Popen(
            [_ffmpeg(), "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
             "-s", "%dx%d" % (w, h), "-r", str(fps), "-i", "-",
             "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", str(crf),
             "-movflags", "+faststart", "-loglevel", "error", out_path],
            stdin=subprocess.PIPE)
        self.out_path = out_path

    # ---------------------------------------------------------------- drawing
    def _chrome(self, d):
        d.rectangle([0, 0, self.W, self.bar_h], fill=PALETTE["bar"])
        for i, c in enumerate(_LIGHTS):
            cx = 24 + i * 23
            d.ellipse([cx - 7, self.bar_h // 2 - 7, cx + 7, self.bar_h // 2 + 7], fill=c)
        tw = self.tfont.getlength(self.title)
        d.text(((self.W - tw) / 2, self.bar_h // 2 - 10), self.title,
               font=self.tfont, fill=PALETTE["gray"])

    def _norm(self, span):
        if len(span) == 3:
            text, ck, bold = span
        else:
            text, ck = span; bold = False
        return text, ck, bold or ck in _BOLD_KEYS

    def _frame(self, active=None, cursor=False):
        im = Image.new("RGB", (self.W, self.H), PALETTE["bg"])
        d = ImageDraw.Draw(im)
        self._chrome(d)
        rows = list(self.lines)
        if active is not None:
            rows = rows + [active]
        rows = rows[-self.max_lines:]
        y = self.top
        last_x = self.pad_x
        for ri, row in enumerate(rows):
            x = self.pad_x
            for span in row:
                text, ck, bold = self._norm(span)
                d.text((x, y), text, font=(self.fontb if bold else self.font),
                       fill=PALETTE.get(ck, PALETTE["fg"]))
                x += self.cw * len(text)
            last_x = x
            y += self.lh
        if cursor:
            cy = self.top + (len(rows) - 1) * self.lh
            d.rectangle([last_x + 1, cy + 2, last_x + 1 + self.cw, cy + self.lh - 4],
                        fill=(205, 208, 205))
        return im

    def _emit(self, active=None, cursor=False):
        self.proc.stdin.write(self._frame(active, cursor).tobytes())

    # ---------------------------------------------------------------- scripting
    def hold(self, secs, cursor=False, active=None):
        for _ in range(max(1, int(round(secs * self.fps)))):
            self._emit(active=active, cursor=cursor)

    def blank(self, n=1, hold=0.0):
        for _ in range(n):
            self.lines.append([])
        if hold:
            self.hold(hold)

    def line(self, spans, hold=0.0):
        """Commit a full line instantly (output)."""
        self.lines.append(list(spans))
        self._emit()
        if hold:
            self.hold(hold)

    def block(self, rows, per=0.06, hold=0.0):
        """Emit several output lines with a small delay between them."""
        for r in rows:
            self.lines.append(list(r))
            self._emit()
            self.hold(per)
        if hold:
            self.hold(hold)

    def scroll(self, rows, per_frame=1, hold=0.0):
        """Fast-scroll a long block (e.g. a file listing)."""
        buf = []
        for i, r in enumerate(rows):
            self.lines.append(list(r))
            if (i + 1) % per_frame == 0:
                self._emit()
        self._emit()
        if hold:
            self.hold(hold)

    def type(self, prefix, typed, cps=40, hold_after=0.5, hold_before=0.25):
        """Type `typed` (list of spans) after the fixed `prefix` spans, char by char."""
        self.hold(hold_before, cursor=True, active=list(prefix))
        chars = [(c, ck, bold) for (t, ck, bold) in map(self._norm, typed) for c in t]
        per = max(1, int(round(cps / self.fps)))
        i = 0
        while i < len(chars):
            i = min(len(chars), i + per)
            # regroup consecutive chars of the same (ck,bold)
            active = list(prefix)
            for c, ck, bold in chars[:i]:
                if active and active[-1][1] == ck and (len(active[-1]) == 3 and active[-1][2] == bold):
                    active[-1] = (active[-1][0] + c, ck, bold)
                else:
                    active.append((c, ck, bold))
            self._emit(active=active, cursor=True)
        self.lines.append(list(prefix) + list(typed))
        self.hold(hold_after, cursor=True)

    def close(self, tail=1.2):
        self.hold(tail)
        self.proc.stdin.close()
        self.proc.wait()
        return self.out_path


# prompt helper --------------------------------------------------------------
def prompt(host="user@oti", path="~/demo"):
    return [(host, "puser"), (":" + path, "ppath"), ("$ ", "psym")]


def still(out_png, title, lines, w=1180, fontsize=20, pad_bottom=18):
    """Render a single terminal still (no video) to a PNG. `lines` is a list of
    span-lists, same format as Termcast."""
    reg, bold = _mono()
    font = ImageFont.truetype(reg, fontsize); fontb = ImageFont.truetype(bold, fontsize)
    tfont = ImageFont.truetype(reg, 17)
    cw = font.getlength("M"); lh = int(round(fontsize * 1.33))
    bar_h = 44; top = bar_h + 16; pad_x = 26
    h = top + len(lines) * lh + pad_bottom
    im = Image.new("RGB", (w, h), PALETTE["bg"]); d = ImageDraw.Draw(im)
    d.rectangle([0, 0, w, bar_h], fill=PALETTE["bar"])
    for i, c in enumerate(_LIGHTS):
        cx = 24 + i * 23
        d.ellipse([cx - 7, bar_h // 2 - 7, cx + 7, bar_h // 2 + 7], fill=c)
    tw = tfont.getlength(title)
    d.text(((w - tw) / 2, bar_h // 2 - 10), title, font=tfont, fill=PALETTE["gray"])
    y = top
    for row in lines:
        x = pad_x
        for span in row:
            text, ck = span[0], span[1]
            b = (len(span) == 3 and span[2]) or ck in _BOLD_KEYS
            d.text((x, y), text, font=(fontb if b else font), fill=PALETTE.get(ck, PALETTE["fg"]))
            x += cw * len(text)
        y += lh
    im.save(out_png)
    return out_png
