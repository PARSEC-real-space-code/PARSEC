#
# Makefile for parsec user's guide
#
# author: Murilo L. Tiago, UTexas, 2005

all: guide.pdf

guide.pdf: guide.tex
	pdflatex $<
	pdflatex $<

clean:
	rm -rf guide.{aux,idx,log,out,toc,synctex.gz}
