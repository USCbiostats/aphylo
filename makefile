news:
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS
check:
	$(MAKE) install && \
		R CMD check --no-stop-on-test-error aphylo*.tar.gz ; rm aphylo*.tar.gz

checkv:
	$(MAKE) install && \
		R CMD check --as-cran --use-valgrind --no-stop-on-test-error aphylo*.tar.gz; rm aphylo*.tar.gz

install:
	cd ..&& R CMD INSTALL --preclean aphylo/&& cd aphylo/

.PHONY: news check checkv install
