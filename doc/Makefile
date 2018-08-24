.PHONY: all clean

all : HeatSolver.pdf

%.pdf : %.md
	pandoc $< -o $@

clean :
	$(RM) *.pdf
