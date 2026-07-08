"""Wolfram Alpha API integration for symbolic math verification.

Primary endpoint: the **LLM API** (`/v1/llm-api`), per Anthropic's own
cookbook (https://platform.claude.com/cookbook/third-party-wolframalpha-using-llm-api).
That endpoint returns plain text optimized for LLM consumption — Wolfram
Language interpretations, related expressions, definite integrals, special
functions, all in one response.

Also exposes the Simple API (`/v2/result`, single short answer) and the
Full API (`/v2/query`, structured JSON pods) for advanced use.

Setup (one-time):
1. Register a free Developer account at https://developer.wolframalpha.com/
2. Create an App at https://developer.wolframalpha.com/access; copy the App ID
3. Set the env var:  export WOLFRAM_APP_ID=YOUR_APP_ID  (e.g. in ~/.zshrc)
4. (Optional) Add WOLFRAM_APP_ID to a .env file in repo root (gitignored).

The free Developer tier provides 2000 queries/month. Each call to `query_llm()`,
`query_plaintext()`, or `query()` costs 1 query.

Usage examples
--------------
>>> from wolfram_alpha import query_llm, query_plaintext, query
>>> # Preferred for our research (rich output, LLM-friendly):
>>> print(query_llm("integrate cos(pi*m*x/2) from (j-1)*L to j*L"))
>>> # One-line numeric answer:
>>> query_plaintext("integrate cos(pi*x/2) from 0 to 1")
>>> # Structured JSON (advanced):
>>> query("inverse symbolic calculator 0.3803027")

Why use this for the Erdős project
----------------------------------
- Closed-form verification of cell-envelope integrals (cos/sin cell bounds).
- Inverse symbolic lookup for μ candidates (cross-check what PSLQ produces).
- Special-function manipulations (e.g., simplifying tail-bound series).
- Definite integrals where sympy times out or produces unreduced output.
"""
from __future__ import annotations

import os
import urllib.parse
import urllib.request
from typing import Optional


_BASE = "https://api.wolframalpha.com/v2"
_LLM_API = "https://www.wolframalpha.com/api/v1/llm-api"


def _app_id() -> str:
    appid = os.environ.get("WOLFRAM_APP_ID")
    if not appid:
        # Fall back to .env file in repo root
        from pathlib import Path
        env_path = Path(__file__).resolve().parents[2] / ".env"
        if env_path.exists():
            for line in env_path.read_text().splitlines():
                if line.startswith("WOLFRAM_APP_ID="):
                    appid = line.split("=", 1)[1].strip().strip('"').strip("'")
                    break
    if not appid:
        raise RuntimeError(
            "WOLFRAM_APP_ID not set. Register at "
            "https://developer.wolframalpha.com/ and either:\n"
            "  export WOLFRAM_APP_ID=YOUR_APP_ID\n"
            "  or add `WOLFRAM_APP_ID=YOUR_APP_ID` to <repo>/.env"
        )
    return appid


def query_llm(input_str: str, timeout: float = 30.0,
              maxchars: int | None = None) -> str:
    """Wolfram Alpha LLM API — purpose-built endpoint for LLM consumption.

    Returns rich plain text including Wolfram Language interpretations,
    related forms, definite values, and contextual notes. This is the
    endpoint recommended by Anthropic's Wolfram Alpha cookbook and is the
    right default for our research workflow.

    `maxchars` (optional, default ~6800): cap the response length.

    Reference: https://products.wolframalpha.com/llm-api/documentation
    """
    appid = _app_id()
    params = {"input": input_str, "appid": appid}
    if maxchars is not None:
        params["maxchars"] = maxchars
    url = f"{_LLM_API}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def query_plaintext(input_str: str, timeout: float = 30.0) -> str:
    """Return Wolfram Alpha's plaintext answer (the 'simple result' API).

    Fast and cheap — uses the /v2/result endpoint which returns a single
    short answer suitable for verification.
    """
    appid = _app_id()
    params = {"appid": appid, "i": input_str, "format": "plaintext"}
    url = f"{_BASE}/result?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def query(input_str: str, timeout: float = 60.0,
          podformat: str = "plaintext") -> dict:
    """Return the full structured Wolfram Alpha response (the /v2/query API).

    Returns a dict with keys like {'pods': [{'title': ..., 'subpods': [...]}]}
    Useful for multi-pod results (e.g., closed form + decimal + alternative
    representations).
    """
    appid = _app_id()
    params = {
        "appid": appid, "input": input_str,
        "format": podformat, "output": "json",
    }
    url = f"{_BASE}/query?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        import json
        return json.loads(resp.read().decode("utf-8"))


def verify_integral(expression: str, var: str, lower: str, upper: str) -> str:
    """Convenience wrapper for definite integral verification.

    Example:
        verify_integral("cos(pi*m*x/2)", "x", "(j-1)*L", "j*L")
    """
    q = f"integrate {expression} d{var} from {lower} to {upper}"
    return query_plaintext(q)


def inverse_symbolic(decimal: str | float, precision: int = 10) -> str:
    """Try to identify a closed form for a decimal number.

    Wolfram Alpha's 'inverse symbolic calculator' / 'identify' feature is
    similar to PSLQ but covers a much wider base library (special constants,
    π combos, gamma values, etc.).

    Example:
        inverse_symbolic("0.3803027")  → may identify it as some
        combination of π, sqrt, log, etc. (if such a clean form exists)
    """
    if isinstance(decimal, float):
        decimal = f"{decimal:.{precision}f}"
    return query_plaintext(f"identify the number {decimal}")


if __name__ == "__main__":
    # Smoke test (requires WOLFRAM_APP_ID set)
    import sys
    try:
        appid = _app_id()
        print(f"App ID detected (len={len(appid)})\n")
        print("=== LLM API (preferred for our workflow) ===\n")
        print("Test 1: integrate cos(πmx/2) on a unit cell — full LLM output:")
        print(query_llm("integrate cos(pi*m*x/2) dx from (j-1)*L to j*L"))
        print("\n" + "=" * 60)
        print("Test 2: inverse-symbolic on our LB headline 0.3803027:")
        print(query_llm("identify the number 0.3803027 closed form"))
        print("\n" + "=" * 60)
        print("Test 3: inverse-symbolic on Together's UB 0.380871:")
        print(query_llm("identify the number 0.380871 closed form"))
    except RuntimeError as e:
        print(f"Setup not complete: {e}")
        sys.exit(1)
