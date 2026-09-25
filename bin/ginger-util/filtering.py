# Copyright (C) 2026 Itoh Laboratory, Institute of Science Tokyo
# 
# This file is part of GINGER.
# 
# GINGER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GINGER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with GINGER; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


import sys
import re
import io


class LookAheadReader:
    def __init__(self, iterable):
        self.iterator = iter(iterable)
        self.peeking = None
        self.count = 0

    def peek(self):
        if self.peeking is None:
            self.peeking = next(self.iterator, None)
        return self.peeking

    def next(self):
        line = self.peek()
        if line is None:
            raise ValueError("unexpected end of input")
        self.peeking = None
        self.count += 1
        return line

    def expect(self, pattern):
        line = self.next()
        m = re.match(pattern, line)
        if not m:
            raise ValueError(f"unexpected pattern at {self.count}: {line!r}")
        return m

    def expect_blank(self):
        line = self.next()
        if line != "\n":
            raise ValueError(f"unexpected non-blank at {self.count}: {line!r}")
        return line


class Path:
    def __init__(self, text, identity):
        self.text = text
        self.identity = identity
        self.alignment = None


class Query:
    def __init__(self, name, pathline):
        self.name = name
        self.pathline = pathline
        self.paths = []

    def write(self, out, threshold):
        append_name = True if len(self.paths) > 1 else False
        for i, path in enumerate(self.paths, start=1):
            if path.identity <= threshold:
                continue
            name = f"{self.name}_path{i}" if append_name else self.name
            out.write(f"{name}\n{self.pathline}{path.text}")
            if path.alignment is not None:
                out.write(f"Alignments:\n{path.alignment}")


def is_alignment_end(line):
    return line is None or line.startswith(("  Alignment for path", ">"))


def read_alignments(reader, paths):
    for i, path in enumerate(paths, start=1):
        lines = [reader.expect(f"  Alignment for path {i}:").string]
        lines.append(reader.expect_blank())
        while not is_alignment_end(reader.peek()):
            lines.append(reader.next())
        path.alignment = "".join(lines)


def read_block(reader):
    lines = []
    while reader.peek() is not None:
        line = reader.next()
        lines.append(line)
        if line == "\n":
            break
    return "".join(lines)


def read_paths(reader, n):
    paths = []
    for i in range(1, n + 1):
        line = reader.peek()
        if line is None or line.startswith(("Alignments:", ">")):
            break
        header = reader.expect(f"  Path {i}:").string
        text = header + read_block(reader)
        m = re.search(r"Percent identity: (\S+) ", text)
        paths.append(Path(text, float(m.group(1)) if m else 0.0))
    return paths


def read_query(reader, name):
    pathline = reader.expect(r"Paths \((\d+)\):")
    n = int(pathline.group(1))
    query = Query(name, pathline.string)
    if n == 0:
        reader.expect_blank()
    query.paths = read_paths(reader, n)
    line = reader.peek()
    if line is not None and line == "Alignments:\n":
        reader.next()
        read_alignments(reader, query.paths)
    return query


def read_queries(reader):
    while reader.peek() is not None:
        m = re.match(r">.+", reader.next())
        if m:
            yield read_query(reader, m.group(0))


def main():
    input_path = sys.argv[1]
    identity_more_than = float(sys.argv[2])
    # "latin-1" to prevent decoding
    out = io.TextIOWrapper(sys.stdout.buffer, encoding="latin-1")
    with open(input_path, encoding="latin-1") as f:
        for query in read_queries(LookAheadReader(f)):
            query.write(out, identity_more_than)
    out.detach()


if __name__ == "__main__":
    main()
