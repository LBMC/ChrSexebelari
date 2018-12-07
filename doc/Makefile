all: TP_experimental_biologists.pdf TP_computational_biologists.pdf

TP_experimental_biologists.pdf: TP_experimental_biologists.md
	R -e 'require(rmarkdown); rmarkdown::render("TP_experimental_biologists.md")'

TP_computational_biologists.pdf: TP_computational_biologists.md
	R -e 'require(rmarkdown); rmarkdown::render("TP_computational_biologists.md")'
