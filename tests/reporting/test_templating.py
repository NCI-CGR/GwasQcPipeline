from dataclasses import dataclass

import pytest
from jinja2.nativetypes import NativeEnvironment


@pytest.fixture(scope="module")
def env():
    return NativeEnvironment()


def test_render(env):
    tpl = env.from_string("{% if test_num > 0 -%}Yes{%- else -%}No{%- endif %}")
    assert "No" == tpl.render(test_num=0)
    assert "Yes" == tpl.render(test_num=5)


def test_render_dataclass(env):
    @dataclass
    class data:
        num: int

    tpl = env.from_string("{% if test.num > 0 -%}Yes{%- else -%}No{%- endif %}")
    assert "No" == tpl.render(test=data(0))
    assert "Yes" == tpl.render(test=data(5))
