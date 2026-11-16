"""Read Shimadzu LabSolutions .lcd data files (Microsoft OLE compound documents), read-only.

What is read, and where it lives inside the file:
  File Property                       sample name/ID/type, vial, tray, injection volume, dilution
  LSS Raw Data/2D Data Item           channel map: which Chromatogram Ch# is which detector model
  LSS Raw Data/Chromatogram Ch#       raw detector signal, delta-encoded integers
  LSS Data Processing Original/PT-*   peak table as processed at acquisition
  LSS Data Processing/PT-*            peak table as it stands now (after any reprocessing)

The layouts were worked out on LabSolutions 5.97 SP1 files and cross-checked against LabSolutions' own exports.
The raw-signal decoding follows chromConverter's `decode_sz_block` (Ethan Bass, GPL-3), re-implemented here.
"""

from __future__ import annotations

import re
import struct
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from pathlib import Path

import numpy as np
import olefile

RID_MODELS = ("RID-20A", "RID-10A", "RID-40")
_PEAK_RECORD = 792  # bytes per peak in a "VER1" peak table
_MANUAL_BIT = 0x8000  # set on peaks drawn or edited by hand (matches the "M" mark in exports)


@dataclass(frozen=True)
class Peak:
    rt: float           # min, apex
    area: float         # µV·s
    height: float       # µV
    start: float        # min
    end: float          # min
    manual: bool


@dataclass(frozen=True)
class Channel:
    number: int         # the N of "Chromatogram ChN"
    source_id: str      # e.g. "LC.1.1.DET.2.CH#1"; the peak table stream is "PT-" + this
    model: str          # e.g. "RID-20A"
    name: str           # e.g. "Detector B Channel 1"
    factor: float       # µV per raw count
    interval_ms: int
    delay_ms: float


@dataclass
class Injection:
    path: str
    sample_name: str
    sample_id: str
    sample_type: str
    vial: str
    tray: str
    inj_vol_uL: float | None
    dil_factor: float | None
    acquired: str               # ISO time, UTC
    modified_count: int
    channels: list[Channel] = field(default_factory=list)
    rid: Channel | None = None
    peaks: list[Peak] = field(default_factory=list)            # current peak table, RID channel
    peaks_original: list[Peak] = field(default_factory=list)   # as processed at acquisition


def _decode_value(v: str):
    """Numbers and strings in LabSolutions XML are hex: @DtoX@ (float64), @FtoX@ (float32), @StoX@ (text).
    Short literals such as "@FtoX@1" stand for the plain value."""
    for tag, size, fmt in (("@DtoX@", 16, ">d"), ("@FtoX@", 8, ">f")):
        if v.startswith(tag):
            h = v[len(tag):]
            if len(h) < 4:
                return float(h or 0)
            return struct.unpack(fmt, bytes.fromhex(h.rjust(size, "0")))[0]
    if v.startswith("@StoX@"):
        return bytes.fromhex(v[6:]).decode("cp1252", "replace")
    return v


def _xml(ole, stream: str) -> str:
    raw = ole.openstream(stream).read()
    return raw[raw.find(b"<"):].decode("cp1252", "replace")


def _tag(xml: str, key: str):
    m = re.search(rf"<{key}>(.*?)</{key}>", xml, re.S)
    return _decode_value(m.group(1)) if m else None


def _filetime(low, high) -> str:
    try:
        ticks = (int(high) << 32) | int(low)
    except (TypeError, ValueError):
        return ""
    if not ticks:
        return ""
    when = datetime(1601, 1, 1, tzinfo=timezone.utc) + timedelta(microseconds=ticks / 10)
    return when.isoformat(timespec="seconds") if when.year >= 2000 else ""  # unset fields decode to ~1601


def _channels(ole) -> list[Channel]:
    if not ole.exists("LSS Raw Data/2D Data Item"):
        return []
    out = []
    for seg in re.findall(r"<DII[^>]*>(.*?)</DII>", _xml(ole, "LSS Raw Data/2D Data Item"), re.S):
        source_id = _tag(seg, "DSID") or ""
        if ".CH#" not in source_id:
            continue  # status traces (pressure, temperatures)

        def le_double(key: str, default: float) -> float:
            m = re.search(rf"<{key}>([0-9A-Fa-f]{{16}})</{key}>", seg)
            return struct.unpack("<d", bytes.fromhex(m.group(1)))[0] if m else default

        out.append(Channel(number=int(_tag(seg, "CN") or 0), source_id=source_id, model=_tag(seg, "DSN") or "",
                           name=_tag(seg, "DN") or "", factor=le_double("CF", 1.0),
                           interval_ms=int(_tag(seg, "Rate") or 500), delay_ms=le_double("DLT", 0.0)))
    return out


def _peak_table(raw: bytes) -> list[Peak]:
    if len(raw) < 20 or raw[:4] != b"VER1":
        return []
    n = struct.unpack_from("<i", raw, 4)[0]
    if n <= 0 or (len(raw) - 20) // n != _PEAK_RECORD:
        return []
    peaks = []
    for k in range(n):
        o = 20 + k * _PEAK_RECORD
        flag, rt_ms = struct.unpack_from("<ii", raw, o)
        area = struct.unpack_from("<d", raw, o + 8)[0]
        height = struct.unpack_from("<d", raw, o + 24)[0]
        start_ms, end_ms = struct.unpack_from("<ii", raw, o + 56)
        peaks.append(Peak(rt=rt_ms / 60000, area=area, height=height, start=start_ms / 60000,
                          end=end_ms / 60000, manual=bool(flag & _MANUAL_BIT)))
    return peaks


def _read(path: str) -> Injection:
    ole = olefile.OleFileIO(path)
    try:
        prop = _xml(ole, "File Property")
        channels = _channels(ole)
        rid = next((c for c in channels if c.model in RID_MODELS), None)
        peaks, original = [], []
        if rid:
            for storage, target in (("LSS Data Processing", peaks), ("LSS Data Processing Original", original)):
                stream = f"{storage}/PT-{rid.source_id}"
                if ole.exists(stream):
                    target.extend(_peak_table(ole.openstream(stream).read()))
        return Injection(
            path=path,
            sample_name=str(_tag(prop, "smpl_name") or ""),
            sample_id=str(_tag(prop, "smpl_id") or ""),
            sample_type=str(_tag(prop, "smpl_type") or ""),
            vial=str(_tag(prop, "vial_num") or ""),
            tray=str(_tag(prop, "tray_name") or ""),
            inj_vol_uL=_tag(prop, "inj_vol"),
            dil_factor=_tag(prop, "dil_factor"),
            acquired=_filetime(_tag(prop, "dwLowDateTime"), _tag(prop, "dwHighDateTime")),
            modified_count=int(_tag(prop, "nModifiedCount") or 0),
            channels=channels, rid=rid, peaks=peaks, peaks_original=original,
        )
    finally:
        ole.close()


@lru_cache(maxsize=4096)
def _read_cached(path: str, mtime: float, size: int) -> Injection:
    return _read(path)


def read_injection(path: str | Path) -> Injection:
    """Parse one .lcd. Cached by (path, mtime, size), so a re-copied file is read again."""
    p = Path(path)
    st = p.stat()
    return _read_cached(str(p), st.st_mtime, st.st_size)


def _deltas(raw: bytes) -> np.ndarray:
    """One delta-encoded sub-block. The high nibble of each record's first byte gives its length and
    sign (even = positive, odd = negative two's complement); 0x82 is a marker, not a value."""
    b = raw + b"\0\0\0\0"
    out, i, n = [], 0, len(raw)
    while i < n:
        x = b[i]
        if x == 0x82:
            i += 1
            continue
        sign_digit = x >> 4
        size = 1 + sign_digit // 2 if sign_digit > 1 else 1
        v = x & 0x0F
        if size == 1:
            v = v - 16 if sign_digit == 1 else v
        else:
            for k in range(1, size):
                v = v * 256 + b[i + k]
            if sign_digit % 2:
                v -= 2 ** (8 * size - 4)
        out.append(v)
        i += size
    return np.cumsum(np.array(out, dtype=np.float64))


def chromatogram(path: str | Path, channel: Channel) -> tuple[np.ndarray, np.ndarray]:
    """Raw signal of one channel: (time in min, signal in µV)."""
    ole = olefile.OleFileIO(str(path))
    try:
        raw = ole.openstream(f"LSS Raw Data/Chromatogram Ch{channel.number}").read()
    finally:
        ole.close()
    interval_ms = struct.unpack_from("<i", raw, 4)[0] or channel.interval_ms
    n_points = struct.unpack_from("<i", raw, 8)[0]
    pos, parts, got = 24, [], 0
    while got < n_points and pos + 2 <= len(raw):
        n_bytes = struct.unpack_from("<h", raw, pos)[0]
        pos += 2
        if n_bytes <= 0:
            break
        block = _deltas(raw[pos:pos + n_bytes])
        parts.append(block)
        got += len(block)
        pos += n_bytes + 2  # skip the trailing length marker
    signal = np.concatenate(parts)[:n_points] if parts else np.zeros(0)
    t = (channel.delay_ms + np.arange(len(signal)) * interval_ms) / 60000
    return t, signal * channel.factor
