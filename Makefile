VERSION:=$(shell Rscript -e 'x<-readLines("DESCRIPTION");cat(gsub(".+[:]\\s*", "", x[grepl("^Vers", x)]))')
PKGNAME:=$(shell Rscript -e 'x<-readLines("DESCRIPTION");cat(gsub(".+[:]\\s*", "", x[grepl("^Package", x)]))')

install: 
	$(MAKE) clean && \
		R CMD build . && \
		R CMD INSTALL $(PKGNAME)_$(VERSION).tar.gz
		

$(PKGNAME)_$(VERSION).tar.gz: R/*.R inst/NEWS README.md
	R CMD build --no-build-vignettes --no-manual . 

inst/NEWS: NEWS.md
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')" && \
	head -n 80 inst/NEWS

README.md: README.qmd
	quarto render README.qmd 

.PHONY: checfull checkv clean

check: $(PKGNAME)_$(VERSION).tar.gz
	R CMD check --no-vignettes --no-manual $(PKGNAME)_$(VERSION).tar.gz

checkfull: R/*.R inst/NEWS README.md
	R CMD build . && \
		R CMD check --as-cran $(PKGNAME)_$(VERSION).tar.gz

checkv: $(PKGNAME)_$(VERSION).tar.gz
	R CMD check --as-cran --use-valgrind $(PKGNAME)_$(VERSION).tar.gz

clean:
	rm -rf $(PKGNAME).Rcheck $(PKGNAME)_$(VERSION).tar.gz; \
		Rscript --vanilla -e 'devtools::clean_dll();devtools::clean_vignettes()'

.PHONY: man docker
man: R/* 
	Rscript --vanilla -e 'roxygen2::roxygenize()'

docker:
	docker run -v$(pwd):/pkg/ -w/pkg --rm -i uscbiostats/fmcmc:latest make check
