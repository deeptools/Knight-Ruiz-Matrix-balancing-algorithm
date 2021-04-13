all : wheels

LABEL:="dilawars/manylinux2010-conan:latest"

build : ./Dockerfile ./build_wheels_linux.sh 
	mkdir -p $(HOME)/wheelhouse
	# Generate version.
	(cd .. && docker build \
	    --rm \
	    -t $(LABEL) \
	    -f .ci/Dockerfile \
	    --build-arg PYPI_PASSWORD=$(PYPI_PASSWORD) . | tee log ) 

wheels : ./Dockerfile ./build_wheels_linux.sh build
	mkdir -p $(HOME)/wheelhouse
	# Generate version.
	docker run -e PYPI_PASSWORD=$(PYPI_PASSWORD) $(LABEL)

test:
	docker run -it $(LABEL) bash
