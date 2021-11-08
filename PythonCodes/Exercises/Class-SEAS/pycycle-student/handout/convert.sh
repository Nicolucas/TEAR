#!/bin/bash
pandoc -f markdown-auto_identifiers -t latex ../README.md -o exercises.tex -s --lua-filter gitlab-math.lua
