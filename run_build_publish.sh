#!/bin/bash


jb build --all --path-output ./book/ ./book/docs

read -p 'Continue (y/n)? ' cont

if [ -z $cont ]
then
	cont='y'
fi



if [ $cont = 'n' ] 
then
	exit
fi

#echo "RUN"

ghp-import -n -p -f ./book/_build/html

#file:///home/krojas/jupyter_book/quantum-espresso-book/book/_build/html/index.html
