"""Vectorised expressions over particle properties.

Accepts Python syntax and C-style operators:
    type == 2 && z > 10
    abs(c_pe) > 3.5 || selection
    sqrt(vx^2 + vy^2 + vz^2)
`and`/`or`/`not` (and && / || / !) work element-wise on arrays. `^` means power.
"""
from __future__ import annotations

import ast
import re

import numpy as np

from .data import Frame, scalar_columns

_FUNCS = {
    name: getattr(np, name)
    for name in ("abs", "sqrt", "exp", "log", "log10", "sin", "cos", "tan", "arcsin", "arccos",
                 "arctan", "arctan2", "floor", "ceil", "round", "minimum", "maximum", "where",
                 "clip", "sign", "mod", "hypot", "isnan", "isfinite")
}
_FUNCS.update({"min": np.minimum, "max": np.maximum, "pi": np.pi, "inf": np.inf, "rint": np.rint,
               "fmod": np.fmod, "np": np})


class _Vectorize(ast.NodeTransformer):
    def visit_BoolOp(self, node):
        self.generic_visit(node)
        fn = "logical_and" if isinstance(node.op, ast.And) else "logical_or"
        out = node.values[0]
        for v in node.values[1:]:
            out = ast.Call(func=ast.Attribute(value=ast.Name("np", ast.Load()), attr=fn, ctx=ast.Load()),
                           args=[out, v], keywords=[])
        return out

    def visit_UnaryOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Not):
            return ast.Call(func=ast.Attribute(value=ast.Name("np", ast.Load()), attr="logical_not", ctx=ast.Load()),
                            args=[node.operand], keywords=[])
        return node

    def visit_Compare(self, node):
        # Chained comparisons (0 < x < 5) -> element-wise and of each comparison.
        self.generic_visit(node)
        if len(node.ops) == 1:
            return node
        parts = []
        left = node.left
        for op, right in zip(node.ops, node.comparators):
            parts.append(ast.Compare(left=left, ops=[op], comparators=[right]))
            left = right
        return self.visit_BoolOp(ast.BoolOp(op=ast.And(), values=parts))

    def visit_BinOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.BitXor):
            node.op = ast.Pow()
        return node

    def visit_Attribute(self, node):
        if not (isinstance(node.value, ast.Name) and node.value.id == "np"):
            raise ValueError("Attribute access is not allowed in expressions")
        return node


def _translate(expr: str) -> str:
    s = expr.replace("&&", " and ").replace("||", " or ")
    s = re.sub(r"!(?!=)", " not ", s)
    # 'Position.X' / 'Color.R' style names -> identifiers
    s = re.sub(r"\b([A-Za-z_]\w*)\.([A-Za-z]\w*)\b",
               lambda m: m.group(0) if m.group(1) == "np" else f"{m.group(1)}_{m.group(2)}", s)
    return s


def compile_expression(expr: str):
    tree = ast.parse(_translate(expr).strip(), mode="eval")
    for node in ast.walk(tree):
        if isinstance(node, (ast.Lambda, ast.ListComp, ast.DictComp, ast.SetComp, ast.GeneratorExp,
                             ast.Subscript, ast.Starred, ast.NamedExpr)):
            raise ValueError(f"'{type(node).__name__}' is not allowed in expressions")
    tree = ast.fix_missing_locations(_Vectorize().visit(tree))
    return compile(tree, "<expression>", "eval")


def namespace(frame: Frame, frame_index: int = 0) -> dict:
    ns = dict(_FUNCS)
    cols = scalar_columns(frame)
    # Case-insensitive aliases so 'Position.X', 'position_x', 'Type' all work.
    for k, v in list(cols.items()):
        ns[k] = v
        ns.setdefault(k.lower(), v)
    for k, v in (("Position_X", "x"), ("Position_Y", "y"), ("Position_Z", "z")):
        ns[k] = ns[k.lower()] = cols[v]
    if "type" in cols:
        ns["ParticleType"] = ns["Particle_Type"] = cols["type"]
    if "id" in cols:
        ns["ParticleIdentifier"] = cols["id"]
    lengths = np.linalg.norm(frame.cell, axis=1)
    ns.update({
        "N": frame.n, "Frame": frame_index, "frame": frame_index,
        "CellSize_X": lengths[0], "CellSize_Y": lengths[1], "CellSize_Z": lengths[2],
        "CellOrigin_X": frame.origin[0], "CellOrigin_Y": frame.origin[1], "CellOrigin_Z": frame.origin[2],
        "ParticleIndex": np.arange(frame.n), "index": np.arange(frame.n),
    })
    for k, v in frame.attributes.items():
        if isinstance(v, (int, float, np.integer, np.floating)):
            ns.setdefault(re.sub(r"\W", "_", k), v)
    return ns


def evaluate(expr: str, frame: Frame, frame_index: int = 0) -> np.ndarray:
    code = compile_expression(expr)
    ns = namespace(frame, frame_index)
    try:
        result = eval(code, {"__builtins__": {}}, ns)
    except NameError as e:
        avail = ", ".join(sorted(scalar_columns(frame)) + ["N", "Frame", "CellSize_X", "…"])
        raise ValueError(f"{e}. Available names: {avail}") from None
    result = np.asarray(result)
    if result.ndim == 0:
        result = np.full(frame.n, result)
    return result


def attribute_namespace(frame: Frame, frame_index: int = 0) -> dict:
    ns = dict(_FUNCS)
    for k, v in frame.attributes.items():
        if isinstance(v, (int, float, np.integer, np.floating)) and not k.startswith("_"):
            ns[re.sub(r"\W", "_", k)] = v
    ns.update(N=frame.n, Frame=frame_index, frame=frame_index)
    ns.setdefault("Timestep", 0)
    return ns


_FIELD = re.compile(r"\{([^{}]+)\}")


def format_overlay(template: str, frame: Frame, frame_index: int = 0) -> str:
    # "{expr}" or "{expr:fmt}" fields evaluated over the frame's numeric attributes, e.g. "t = {Timestep*0.002:.1f} ps".
    ns = attribute_namespace(frame, frame_index)

    def field(m):
        expr, _, spec = _split_spec(m.group(1))
        try:
            v = eval(compile_expression(expr), {"__builtins__": {}}, ns)
            return format(v.item() if hasattr(v, "item") else v, spec)
        except Exception as e:
            return f"<{expr.strip()}: {type(e).__name__}>"
    return _FIELD.sub(field, template)


def _split_spec(s):
    # The format spec follows the last ':' when it looks like one (".3f", "8d", ",") rather than part of the expression.
    i = s.rfind(":")
    if i > 0 and re.fullmatch(r"[<>^=+\- #0-9,._a-zA-Z%]*", s[i + 1:]) and not re.search(r"[A-Za-z_]\w*\(", s[i + 1:]):
        return s[:i], ":", s[i + 1:]
    return s, "", ""
