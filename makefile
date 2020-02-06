

news:
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS

check: aphylo.tar.gz
	R CMD check --no-stop-on-test-error aphylo.tar.gz &

checkv: aphylo.tar.gz
	R CMD check --as-cran --use-valgrind --no-stop-on-test-error aphylo.tar.gz &

aphylo.tar.gz: R/*.R
	rm aphylo.tar.gz; \
	R CMD build . && mv aphylo*.tar.gz aphylo.tar.gz

install: aphylo.tar.gz
	R CMD INSTALL --preclean aphylo.tar.gz

.PHONY: news check checkv install
