news:
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS
check:
	cd ..&&R CMD build aphylo/ && \
		R CMD check --no-stop-on-test-error aphylo*.tar.gz ; rm aphylo*.tar.gz

checkv:
	cd ..&&R CMD build aphylo/ && \
		R CMD check --as-cran --use-valgrind --no-stop-on-test-error aphylo*.tar.gz; rm aphylo*.tar.gz
