#!/bin/bash


jb build --all --path-output ./book/ ./book/docs

ghp-import -n -p -f ./book/_build/html

#file:///home/krojas/jupyter_book/quantum-espresso-book/book/_build/html/index.html