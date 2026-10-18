#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A fluid whose molecules have no partial charges is skipped cleanly."""

from contextlib import contextmanager
from types import SimpleNamespace

import numpy as np
import pytest

import dielectric_constant_step.dielectric_constant as module


class _Printer:
    def __init__(self):
        self.lines = []

    def normal(self, text):
        self.lines.append(str(text))


def _setup(monkeypatch, charges):
    md = SimpleNamespace(charges=lambda group: charges)

    @contextmanager
    def open_md(path):
        yield md

    monkeypatch.setattr(module, "open_md", open_md)
    printer = _Printer()
    monkeypatch.setattr(module, "printer", printer)
    n = len(charges)
    configuration = SimpleNamespace(
        find_molecules=lambda as_indices: [[i] for i in range(n)],
        atoms=SimpleNamespace(symbols=["Ar"] * n),
    )
    step = SimpleNamespace(paths=["run1"], title="Dielectric Constant", indent="")
    return step, configuration, printer


def test_argon_is_skipped(monkeypatch):
    step, configuration, printer = _setup(monkeypatch, [0.0] * 4)
    result = module.DielectricConstant.analyze(
        step, P={"dipole of": "automatic"}, configuration=configuration
    )
    assert result is None
    text = " ".join(printer.lines)
    assert "no partial charges" in text and "skipped" in text


def test_charged_molecules_are_not_skipped(monkeypatch):
    """Atoms with charges go on to the analysis (here, to reading the runs)."""
    step, configuration, printer = _setup(monkeypatch, [0.4, -0.4, 0.4, -0.4])
    configuration.find_molecules = lambda as_indices: [[0, 1], [2, 3]]
    step._read = lambda *args: (_ for _ in ()).throw(RuntimeError("reading"))
    with pytest.raises(RuntimeError, match="reading"):
        module.DielectricConstant.analyze(
            step, P={"dipole of": "automatic"}, configuration=configuration
        )
    assert not any("no partial charges" in line for line in printer.lines)


def test_no_charge_threshold():
    assert np.all(np.abs(np.zeros(3)) < module.NO_CHARGE)
