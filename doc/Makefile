all: pdf

R_OPTS=--no-save --no-restore --no-init-file --no-site-file

pdf: pred_wt.R
	cp pred_wt.R pred_wt.Rmd
	R ${R_OPTS} -e 'library(knitr); knit("pred_wt.Rmd"); knit("app.Rmd")'
	pandoc -V geometry:margin=1in pred_wt.md -o pred_wt.pdf --bibliography=../pred_wt.bib
	pandoc -V geometry:margin=1in app.md -o app.pdf 
	pandoc -V geometry:margin=1in pred_wt.md app.md -o pred_wt.html --bibliography=../references/Gmacs.bib

clean:
	rm pred_wt.md 
	rm app.md 
