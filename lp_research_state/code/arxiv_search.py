"""arXiv API search — free, no auth, callable directly via stdlib.

Use this from Bash (or import as a module) when you need to find papers
on minimum overlap, autocorrelation bounds, or related combinatorics work
that has appeared since our last literature scan.

The arXiv export API is documented at https://arxiv.org/help/api/user-manual.
No key required; rate limit is 1 request per 3 seconds (be polite).
"""
from __future__ import annotations

import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Iterator


_API = "http://export.arxiv.org/api/query"
_NS = {"atom": "http://www.w3.org/2005/Atom",
       "arxiv": "http://arxiv.org/schemas/atom"}


@dataclass
class Paper:
    arxiv_id: str
    title: str
    authors: list[str]
    abstract: str
    published: str   # YYYY-MM-DD
    categories: list[str]
    pdf_url: str
    abs_url: str

    def short(self) -> str:
        a = ", ".join(self.authors[:3])
        if len(self.authors) > 3:
            a += f", +{len(self.authors)-3} more"
        return f"[{self.arxiv_id}] {self.published}\n  {self.title}\n  {a}"


def search(query: str, max_results: int = 20,
           sort_by: str = "submittedDate",
           sort_order: str = "descending",
           start: int = 0) -> Iterator[Paper]:
    """Run an arXiv search and yield Paper objects."""
    params = {
        "search_query": query,
        "start": start,
        "max_results": max_results,
        "sortBy": sort_by,
        "sortOrder": sort_order,
    }
    url = f"{_API}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = resp.read()
    root = ET.fromstring(data)
    for entry in root.findall("atom:entry", _NS):
        link_pdf = ""
        link_abs = ""
        for link in entry.findall("atom:link", _NS):
            if link.attrib.get("title") == "pdf":
                link_pdf = link.attrib.get("href", "")
            elif link.attrib.get("rel") == "alternate":
                link_abs = link.attrib.get("href", "")
        arxiv_id = ""
        id_el = entry.find("atom:id", _NS)
        if id_el is not None and id_el.text:
            arxiv_id = id_el.text.rsplit("/", 1)[-1]
        title = (entry.find("atom:title", _NS).text or "").strip()
        abstract = (entry.find("atom:summary", _NS).text or "").strip()
        published = (entry.find("atom:published", _NS).text or "")[:10]
        authors = [
            (a.find("atom:name", _NS).text or "").strip()
            for a in entry.findall("atom:author", _NS)
        ]
        cats = [c.attrib.get("term", "")
                for c in entry.findall("atom:category", _NS)]
        yield Paper(arxiv_id=arxiv_id, title=title, authors=authors,
                    abstract=abstract, published=published,
                    categories=cats, pdf_url=link_pdf, abs_url=link_abs)


def search_min_overlap(years: tuple[int, int] = (2024, 2026),
                       max_results: int = 20) -> list[Paper]:
    """Search for recent papers on Erdős minimum overlap or related topics."""
    q = (
        'abs:"minimum overlap" OR abs:"Erdős overlap" '
        'OR abs:"autocorrelation lower bound" OR ti:"minimum overlap"'
    )
    results = list(search(q, max_results=max_results))
    y0, y1 = years
    return [p for p in results
            if y0 <= int(p.published[:4]) <= y1]


if __name__ == "__main__":
    import sys
    query = sys.argv[1] if len(sys.argv) > 1 else 'abs:"minimum overlap" AND abs:"Erdős"'
    print(f"Searching arXiv: {query}\n")
    n = 0
    for p in search(query, max_results=15):
        print(p.short())
        print()
        n += 1
    print(f"({n} results)")
