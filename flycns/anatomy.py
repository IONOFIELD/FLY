"""Braille dorsal-view projection of the CNS and frame rendering, shared by viz tools."""
import os, sys, time
import numpy as np
import pandas as pd

PAL = {"ol_intrinsic": 36, "ol_sensory": 36, "visual_projection": 96, "visual_centrifugal": 96,
       "cb_intrinsic": 33, "cb_sensory": 93, "cb_motor": 91, "descending_neuron": 35,
       "ascending_neuron": 95, "vnc_intrinsic": 32, "vnc_motor": 31, "vnc_sensory": 92, "vnc_efferent": 91}
MARK = {"DNp01": "G", "TTMn": "T", "PSI": "P", "MN9": "M"}
BR = [[0x01, 0x08], [0x02, 0x10], [0x04, 0x20], [0x40, 0x80]]
DIM, RESET = "\033[2;37m", "\033[0m"


class Projector:
    def __init__(self, neurons, W=None, H=None, color=None, landscape=None, stretch=False):
        tty = sys.stdout.isatty()
        if tty:
            ts = os.get_terminal_size()
        self.W = W or (min(ts.columns - 1, 240) if tty else 120)
        self.H = H or (max(ts.lines - 6, 20) if tty else 40)
        # landscape: brain on the left, VNC to the right (fills wide terminals)
        self.landscape = (self.W * 2 > self.H * 4 * 1.2) if landscape is None else landscape
        self.stretch = stretch
        self.color = tty if color is None else color
        self.tty = tty
        n = neurons.dropna(subset=["x", "y", "z"]).copy()
        X = n[["x", "y", "z"]].values.astype(float)
        sc = n["superclass"].fillna("")
        brain = X[sc.str.startswith(("cb_", "ol_")).values].mean(0)
        vnc = X[sc.str.startswith("vnc").values].mean(0)
        ap = vnc - brain; ap /= np.linalg.norm(ap)
        side = n["somaSide"].values
        lr = X[side == "R"].mean(0) - X[side == "L"].mean(0)
        lr -= ap * (lr @ ap); lr /= np.linalg.norm(lr)
        self.brain, self.ap, self.lr = brain, ap, lr
        self.DW, self.DH = self.W * 2, self.H * 4
        u = (X - brain) @ lr; v = (X - brain) @ ap
        if self.landscape:          # horizontal = brain->VNC, vertical = animal's right->left
            u, v = v, -u
        self.ulo, uhi = np.percentile(u, [0.3, 99.7]); self.vlo, vhi = np.percentile(v, [0.3, 99.7])
        if self.stretch:
            self.sx, self.sy = self.DW / (uhi - self.ulo), self.DH / (vhi - self.vlo)
        else:
            s0 = min(self.DW / (uhi - self.ulo), self.DH / (vhi - self.vlo)); self.sx = self.sy = s0
        self.ox = (self.DW - int((uhi - self.ulo) * self.sx)) // 2
        self.oy = (self.DH - int((vhi - self.vlo) * self.sy)) // 2
        px, py = self.project(X)
        n["dx"], n["dy"] = px, py
        self.pos = n.set_index("bodyId")[["dx", "dy"]]
        self.meta = n.set_index("bodyId")
        dens = np.zeros((self.DH, self.DW), int); np.add.at(dens, (py, px), 1)
        thr = np.quantile(dens[dens > 0], 0.35)
        self.bg = dens >= thr

    def project(self, xyz):
        u = (xyz - self.brain) @ self.lr; v = (xyz - self.brain) @ self.ap
        if self.landscape:
            u, v = v, -u
        px = np.clip(((u - self.ulo) * self.sx).astype(int) + self.ox, 0, self.DW - 1)
        py = np.clip(((v - self.vlo) * self.sy).astype(int) + self.oy, 0, self.DH - 1)
        return px, py

    @property
    def view_label(self):
        return "brain LEFT, VNC RIGHT" if self.landscape else "brain TOP, VNC BOTTOM"

    def cell(self, mask, cy, cx):
        bits = 0
        for r in range(4):
            for c in range(2):
                if mask[cy * 4 + r, cx * 2 + c]:
                    bits |= BR[r][c]
        return chr(0x2800 + bits) if bits else " "

    def col(self, code, s, bold=False):
        return f"\033[{'1;' if bold else ''}{code}m{s}{RESET}" if self.color else s

    def base_grid(self, layer=None, outline=True):
        layer = self.bg if layer is None else layer
        g = [[" "] * self.W for _ in range(self.H)]
        for cy in range(self.H):
            for cx in range(self.W):
                ch = self.cell(layer, cy, cx)
                if ch != " ":
                    g[cy][cx] = (DIM + ch + RESET) if self.color else ch
                elif outline and layer is not self.bg:
                    ch2 = self.cell(self.bg, cy, cx)
                    if ch2 != " ":
                        g[cy][cx] = (DIM + ch2 + RESET) if self.color else "."
        return g

    def animate(self, spikes, onset, window_ms=200, bin_ms=5, fade=3, title="", base=None, links=None, delay=0.2):
        """spikes: DataFrame with t_ms, bodyId, type. Plays frames from onset."""
        w = spikes[(spikes.t_ms >= onset) & (spikes.t_ms < onset + window_ms) & spikes.bodyId.isin(self.pos.index)].copy()
        w["b"] = ((w.t_ms - onset) // bin_ms).astype(int)
        base = base or self.base_grid()
        for b in range(window_ms // bin_ms):
            grid = [r[:] for r in base]
            act = np.zeros((self.DH, self.DW), bool); actsc = {}; letters = {}
            for age in range(fade, -1, -1):
                cur = w[w.b == b - age]
                for bid, t in zip(cur.bodyId, cur.type):
                    x, y = self.pos.at[bid, "dx"], self.pos.at[bid, "dy"]
                    s = self.meta.at[bid, "superclass"]
                    if t in MARK:
                        letters[(y // 4, x // 2)] = (MARK[t], s)
                    else:
                        act[y, x] = True; actsc[(y // 4, x // 2)] = (s, age)
            if links is not None:
                firing = set(w[w.b == b].bodyId)
                if firing:
                    L = links[links.bodyId_pre.isin(firing)]
                    syn = np.zeros((self.DH, self.DW), bool); syn[L.dy.values, L.dx.values] = True
                    for cy in range(self.H):
                        for cx in range(self.W):
                            ch = self.cell(syn, cy, cx)
                            if ch != " ":
                                grid[cy][cx] = self.col(93, ch, True)
            for (cy, cx), (s, age) in actsc.items():
                ch = self.cell(act, cy, cx)
                if ch != " ":
                    grid[cy][cx] = self.col(PAL.get(s, 37), ch, bold=(age == 0))
            for (cy, cx), (ch, s) in letters.items():
                grid[cy][cx] = self.col(PAL.get(s, 37), ch, True)
            if self.tty:
                print("\033[H\033[J", end="")
            print(f" {title}  t={b*bin_ms:>3d} ms  {self.view_label}  [G GF  T TTMn  P PSI  M MN9]")
            for r in grid:
                print("".join(r))
            sys.stdout.flush()
            time.sleep(delay if self.tty else 0)
