news:
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS
check:
	$(MAKE) build && \
		R CMD check --no-stop-on-test-error aphylo*.tar.gz ; rm aphylo*.tar.gz

checkv:
	$(MAKE) build && cd ../ && \
		R CMD check --as-cran --use-valgrind --no-stop-on-test-error aphylo*.tar.gz; rm aphylo*.tar.gz 

build:
	cd ../ && R CMD build aphylo/ && cd aphylo/

install:
	cd ..&& R CMD INSTALL --preclean aphylo/&& cd aphylo/

pruner:
	cp ../pruner/include/* inst/include

.PHONY: news check checkv install
