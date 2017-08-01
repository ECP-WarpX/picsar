SRC_FILES = theory.tex \
            Deposition/Current_deposition.tex \
            Gather/Field_gather.tex \
            Maxwell_solvers/Maxwell_FDTD_solver.tex \
            Maxwell_solvers/Maxwell_NSFDTD_solver.tex \
            Maxwell_solvers/Maxwell_PSATD_solver.tex

all: $(SRC_FILES)
	rm -f allbibs.bib;cat *.bib > allbibs.bib
	rm -f theory.md
	pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o theorycore.md
	python correct.py