#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `dielectric_constant_step` package."""

import pytest  # noqa: F401
import dielectric_constant_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = dielectric_constant_step.DielectricConstant()
    assert (
        str(type(result))
        == "<class 'dielectric_constant_step.dielectric_constant.DielectricConstant'>"
    )
