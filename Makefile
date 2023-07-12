include .entangled/makefile.in

.PHONY: figures

figures: figures.mk
> make -f figures.mk

