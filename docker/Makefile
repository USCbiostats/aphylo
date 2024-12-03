ifndef CNTR
	CNTR = docker
endif

build:
	$(CNTR) build -t uscbiostats/aphylo:release -f Dockerfile . && \
		$(CNTR) build -t uscbiostats/aphylo:devel -f Dockerfile.dev .

clang:
	$(CNTR) build -t uscbiostats/aphylo:clang -f Dockerfile-clang .

flexiblas:
	$(CNTR) build -t uscbiostats/aphylo:flexiblas -f Dockerfile.flexiblas .

push: build
	$(CNTR) push uscbiostats/aphylo:release && \
		$(CNTR) push uscbiostats/aphylo:devel
