build:
	docker build -t uscbiostats/aphylo:release -f Dockerfile . && \
		docker build -t uscbiostats/aphylo:devel -f Dockerfile.dev .

clang:
	docker build -t uscbiostats/aphylo:clang -f Dockerfile-clang .
	
push: build
	docker push uscbiostats/aphylo:release && \
		docker push uscbiostats/aphylo:devel
