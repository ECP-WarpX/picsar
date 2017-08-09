SRC_FILES = theory.tex \
            Deposition/Current_deposition.tex \
            Gather/Field_gather.tex \
            Maxwell_solvers/Maxwell_FDTD_solver.tex \
            Maxwell_solvers/Maxwell_NSFDTD_solver.tex \
            Maxwell_solvers/Maxwell_PSATD_solver.tex

all: $(SRC_FILES) clean
	pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o theory.md
	python reformat.py

clean: 
	rm -f theory.md
	rm -f theory.txt