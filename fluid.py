#!/usr/bin/env python3
"""
fluid-cinema — JSON tanımından 2B Navier-Stokes akışkan simülasyonu.
Mürekkep, alev, girdap, rüzgar tüneli... hepsi parametre.

Çekirdek: Jos Stam "Stable Fluids" (yarım-Lagrange adveksiyon + basınç
projeksiyonu) + girdap yakalama (vorticity confinement) + yüzdürme.
Çıktı: kare kare PNG -> ffmpeg ile GIF/MP4.

Gereksinim: numpy (render için yalnız stdlib + ffmpeg).
    python3 fluid.py --scene scenes/ink.json
    python3 fluid.py --scene scenes/ink.json --preview
"""

import argparse
import json
import math
import os
import shutil
import struct
import subprocess
import zlib

import numpy as np


# ---------------------------------------------------------------- png

def png_yaz(path, img):
    """img: (H,W,3) uint8 -> RGB PNG (color type 2), yalnız stdlib."""
    h, w, _ = img.shape
    raw = bytearray()
    for y in range(h):
        raw.append(0)
        raw += img[y].tobytes()
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr)
                + chunk(b"IDAT", zlib.compress(bytes(raw), 6))
                + chunk(b"IEND", b""))


# ---------------------------------------------------------------- çözücü

class Akiskan:
    """Stable Fluids + vorticity confinement + buoyancy."""

    def __init__(self, rows, cols, cfg):
        self.rows, self.cols = rows, cols
        self.u = np.zeros((rows, cols))
        self.v = np.zeros((rows, cols))
        self.dye = np.zeros((rows, cols, 3))          # RGB boya
        self.temp = np.zeros((rows, cols))            # yüzdürme sıcaklığı
        self.solid = np.zeros((rows, cols), bool)     # engeller
        self.viscosity = cfg.get("viscosity", 0.0005)
        self.dye_diff = cfg.get("dye_diffusion", 0.00001)
        self.dissipation = cfg.get("dissipation", 0.996)
        self.vorticity = cfg.get("vorticity", 12.0)
        self.buoyancy = cfg.get("buoyancy", 0.0)
        self.pressure_iters = cfg.get("pressure_iters", 40)
        self.smoothing = cfg.get("velocity_smoothing", 0.12)

    # ---- yardımcılar
    def _diffuse(self, f, k):
        if k <= 0:
            return f
        lap = np.zeros_like(f)
        lap[1:-1, 1:-1] = (f[:-2, 1:-1] + f[2:, 1:-1] + f[1:-1, :-2]
                           + f[1:-1, 2:] - 4 * f[1:-1, 1:-1])
        return f + k * lap

    def _advect(self, f, dt):
        h, w = f.shape[:2]
        yy, xx = np.mgrid[0:h, 0:w].astype(np.float64)
        bx = self._yansit(xx - self.u * dt, w)
        by = self._yansit(yy - self.v * dt, h)
        x0 = bx.astype(int); y0 = by.astype(int)
        x1 = np.minimum(x0 + 1, w - 1); y1 = np.minimum(y0 + 1, h - 1)
        sx = bx - x0; sy = by - y0
        if f.ndim == 3:
            sx = sx[..., None]; sy = sy[..., None]
        return (f[y0, x0] * (1 - sx) * (1 - sy) + f[y0, x1] * sx * (1 - sy)
                + f[y1, x0] * (1 - sx) * sy + f[y1, x1] * sx * sy)

    @staticmethod
    def _yansit(t, n):
        """CFL güvenli geri-izleme: örnekleyici [1, n-2] içinde kalır."""
        return np.clip(t, 1.0, n - 2.001)

    def _project(self):
        h, w = self.u.shape
        div = np.zeros((h, w))
        p = np.zeros((h, w))
        div[1:-1, 1:-1] = -0.5 * (self.u[2:, 1:-1] - self.u[:-2, 1:-1]
                                  + self.v[1:-1, 2:] - self.v[1:-1, :-2])
        for _ in range(self.pressure_iters):
            p[1:-1, 1:-1] = (p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2]
                             + p[1:-1, 2:] + div[1:-1, 1:-1]) / 4.0
        self.u[1:-1, 1:-1] -= 0.5 * (p[2:, 1:-1] - p[:-2, 1:-1])
        self.v[1:-1, 1:-1] -= 0.5 * (p[1:-1, 2:] - p[1:-1, :-2])

    def _vorticity(self, dt):
        e = self.vorticity
        if e <= 0:
            return
        h, w = self.u.shape
        wcurl = np.zeros((h, w))
        wcurl[1:-1, 1:-1] = 0.5 * ((self.v[1:-1, 2:] - self.v[1:-1, :-2])
                                   - (self.u[2:, 1:-1] - self.u[:-2, 1:-1]))
        gw_y, gw_x = np.gradient(np.abs(wcurl))
        n = np.sqrt(gw_x**2 + gw_y**2) + 1e-5
        fx = e * (gw_y / n) * wcurl
        fy = -e * (gw_x / n) * wcurl
        self.u += dt * fx
        self.v += dt * fy

    def step(self, dt=1.0):
        # kuvvetler
        self.u = self._diffuse(self.u, self.viscosity)
        self.v = self._diffuse(self.v, self.viscosity)
        if self.buoyancy > 0:
            self.v -= self.buoyancy * self.temp
        self._vorticity(dt)
        np.clip(self.u, -8, 8, out=self.u)
        np.clip(self.v, -8, 8, out=self.v)

        # CFL alt-adımları: en hızlı hücre adım başına ~1.5 hücre yürüsün
        maxv = float(max(np.abs(self.u).max(), np.abs(self.v).max()))
        n = max(1, math.ceil(maxv / 1.5))
        sdt = dt / n
        for _ in range(n):
            self._project()
            self.u = self._advect(self.u, sdt)
            self.v = self._advect(self.v, sdt)
            self.dye = self._advect(self.dye, sdt)
            self.temp = self._advect(self.temp, sdt)
            self._project()
            self._duzelt(self.smoothing)
            # kapalı kutu duvarları
            self.u[:, :3] = self.u[:, -3:] = 0.0
            self.v[:3, :] = self.v[-3:, :] = 0.0

        self.dye *= self.dissipation
        self.temp *= self.dissipation
        for f in (self.u, self.v):
            f[self.solid] = 0.0
        self.dye[self.solid] = 0.0

    def _duzelt(self, k):
        """Damalı kararsızlığı bastıran hafif Laplacian düzeltmesi."""
        if k <= 0:
            return
        for f in (self.u, self.v):
            lap = np.zeros_like(f)
            lap[1:-1, 1:-1] = (f[:-2, 1:-1] + f[2:, 1:-1] + f[1:-1, :-2]
                               + f[1:-1, 2:] - 4.0 * f[1:-1, 1:-1])
            f += k * lap


# ---------------------------------------------------------------- sahne

def hucre(cfg_norm, rows, cols):
    return int(cfg_norm[1] * rows), int(cfg_norm[0] * cols)


def disksa_maske(rows, cols, cy, cx, r):
    yy, xx = np.mgrid[0:rows, 0:cols]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r ** 2


def sahne_uygula(sim, cfg, frame, steps):
    sc = cfg.get("sources", [])
    rows, cols = sim.rows, sim.cols
    for s in sc:
        cy, cx = hucre(s["pos"], rows, cols)
        r = max(1, int(s.get("radius", 0.03) * min(rows, cols)))
        maske = disksa_maske(rows, cols, cy, cx, r)
        if sim.solid[maske].any():
            continue
        vel = s.get("velocity", [0, 0])
        renk = np.array(s.get("color", [1, 1, 1]))
        guc = s.get("strength", 1.0)
        tur = s.get("emit", "jet")
        if tur == "burst" and frame > 0:
            continue
        # hızlar hücre/adım biriminde (CFL ~ 1 civarı)
        sim.u[maske] += vel[0] * guc
        sim.v[maske] += vel[1] * guc
        miktar = s.get("rate", 1.0)
        sim.dye[maske] += renk * miktar * (0.05 if tur == "smolder" else (1.5 if tur == "burst" else 0.18))
        if s.get("temperature", 0) > 0:
            sim.temp[maske] += s["temperature"] * 0.3

    for f in cfg.get("forces", []):
        t = f["type"]
        if t == "vortex":
            cy, cx = hucre(f["pos"], rows, cols)
            r = max(2, int(f.get("radius", 0.2) * min(rows, cols)))
            maske = disksa_maske(rows, cols, cy, cx, r)
            yy, xx = np.mgrid[0:rows, 0:cols]
            dx, dy = xx - cx, yy - cy
            uzak = np.sqrt(dx**2 + dy**2) + 1e-5
            kacak = f.get("strength", 0.4)
            sim.u[maske] += -dy[maske] / uzak[maske] * kacak
            sim.v[maske] += dx[maske] / uzak[maske] * kacak
        elif t == "wind":
            guc = f.get("strength", 1.0) * 0.1
            bolge = f.get("region")  # [y0,x0,y1,x1] normalize; yoksa hepsi
            if bolge:
                alt = hucre(bolge[:2], rows, cols)
                ust = hucre(bolge[2:], rows, cols)
                sim.u[alt[0]:ust[0], alt[1]:ust[1]] += f["dir"][0] * guc
                sim.v[alt[0]:ust[0], alt[1]:ust[1]] += f["dir"][1] * guc
            else:
                sim.u += f["dir"][0] * guc
                sim.v += f["dir"][1] * guc

    for o in cfg.get("obstacles", []):
        cy, cx = hucre(o["pos"], rows, cols)
        r = max(1, int(o.get("radius", 0.05) * min(rows, cols)))
        sim.solid |= disksa_maske(rows, cols, cy, cx, r)


# ---------------------------------------------------------------- render

def kare_ciz(sim, bg, exposure, vignette, upscale):
    dye = sim.dye
    img = 1.0 - np.exp(-exposure * np.maximum(dye, 0.0))
    img = np.clip(img, 0, 1) ** (1 / 2.2)                      # gamma
    renk = img * 255
    if upscale > 1:
        renk = renk.repeat(upscale, 0).repeat(upscale, 1)
        # 5-nokta kutu bulanıklığı (enerji koruyan, hayaletsiz)
        renk = (renk
                + np.roll(renk, 1, 0) + np.roll(renk, -1, 0)
                + np.roll(renk, 1, 1) + np.roll(renk, -1, 1)) / 5.0
    arka = np.zeros_like(renk)
    arka[..., 0] = bg[0]; arka[..., 1] = bg[1]; arka[..., 2] = bg[2]
    h, w = renk.shape[:2]
    if vignette > 0:
        yy, xx = np.mgrid[0:h, 0:w]
        d = np.sqrt(((yy / h - 0.5) * 2) ** 2 + ((xx / w - 0.5) * 2) ** 2)
        katsayi = 1.0 - vignette * np.clip(d - 0.55, 0, 1) ** 2
        arka *= katsayi[..., None]
    karisim = 1.0 - np.clip(renk.mean(axis=2) / 255.0, 0, 1)
    out = renk + arka * karisim[..., None] * 0.35
    return np.clip(out, 0, 255).astype(np.uint8)


# ---------------------------------------------------------------- ana

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene", default="scenes/ink.json")
    ap.add_argument("--preview", action="store_true")
    args = ap.parse_args()

    cfg = json.load(open(args.scene, encoding="utf-8"))
    sim_cfg = dict(cfg["sim"])
    if args.preview:
        sim_cfg["cols"] = max(48, sim_cfg["cols"] // 3)
        sim_cfg["rows"] = max(27, sim_cfg["rows"] // 3)
        sim_cfg["steps"] = min(sim_cfg["steps"], 80)
        sim_cfg["upscale"] = 2

    cols, rows = sim_cfg["cols"], sim_cfg["rows"]
    steps = sim_cfg["steps"]
    upscale = sim_cfg.get("upscale", 4)
    fps = sim_cfg.get("fps", 30)

    sim = Akiskan(rows, cols, sim_cfg)
    sahne_uygula(sim, cfg, 0, steps)  # statik engelleri kur
    out = cfg["output"]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    tmp = out + ".frames"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)

    bg = cfg.get("background", {})
    bg_renk = bg.get("color", [8, 10, 16])
    exposure = bg.get("exposure", 1.6)
    vignette = bg.get("vignette", 0.4)

    for frame in range(steps):
        sahne_uygula(sim, cfg, frame, steps)
        sim.step(dt=1.0)
        img = kare_ciz(sim, bg_renk, exposure, vignette, upscale)
        png_yaz(f"{tmp}/f{frame:05d}.png", img)
        if frame % 40 == 0:
            print(f"  kare {frame}/{steps}")

    shutil.copy(f"{tmp}/f{steps - 1:05d}.png", out.replace(".gif", ".png"))

    try:
        subprocess.run(
            ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
             "-i", f"{tmp}/f%05d.png",
             "-vf", "palettegen=max_colors=256:stats_mode=diff", f"{tmp}/pal.png"],
            check=True)
        subprocess.run(
            ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
             "-i", f"{tmp}/f%05d.png", "-i", f"{tmp}/pal.png",
             "-lavfi", "paletteuse=dither=bayer:bayer_scale=3:diff_mode=rectangle",
             "-loop", "0", out],
            check=True)
        print(f"== BİTTİ -> {out}")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        print(f"  [uyarı] ffmpeg: {e}; kareler {tmp}/ içinde")
        return
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
