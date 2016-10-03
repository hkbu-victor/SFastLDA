The files included in this directory are intended to help JMLR authors
produce conforming documents using the LyX document preparation
system. Lyx (http://www.lyx.org) is a free system that acts as a
GUI-based front end for LaTeX.

Note that the layout file assumes that both the "ae" and "theapa"
packages are available in the local LaTeX installation.

This directory contains the following files:

README: this file.

jmlr.layout:
algorithmics.inc: layout files which should be placed in ~/.lyx/layouts

cua.bind: emacs-like keybindings for LyX, which should be placed in
	~/.lyx/bind

jmlr2e-lyx.sty: a slight variation on the standard jmlr2e.sty in which
	the initial epsfig, amssymb, natbib and graphicx includes have
	been commented out for LyX compatibility.

jmlr-sample.tgz: an archived sample of a paper prepared using LyX.

These files were graciously prepared and provided by Ralf Herbrich, of
Microsoft Research.
