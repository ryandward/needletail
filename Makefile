IMAGE   ?= needletail-server
TAG     ?= latest
PORT    ?= 8002

.PHONY: docker docker-run docker-push docker-prune build test clean

## Build the Docker image, pruning dangling images left by previous builds
docker:
	docker build -t $(IMAGE):$(TAG) .
	@echo "── pruning dangling images ──"
	docker image prune -f

## Run the container (interactive, auto-remove, detached)
docker-run: docker
	@# Stop any existing container on the same name
	-docker rm -f $(IMAGE) 2>/dev/null
	docker run --rm -d --name $(IMAGE) -p $(PORT):8002 $(IMAGE):$(TAG)
	@echo "$(IMAGE) running on :$(PORT)  (docker logs -f $(IMAGE))"

## Stop the running container
docker-stop:
	-docker rm -f $(IMAGE) 2>/dev/null

## Push to a registry (set IMAGE=registry/repo first)
docker-push: docker
	docker push $(IMAGE):$(TAG)

## Remove all needletail images, build cache, and dangling layers
docker-prune:
	-docker rm -f $(IMAGE) 2>/dev/null
	-docker rmi -f $$(docker images $(IMAGE) -q) 2>/dev/null
	docker builder prune -f
	docker image prune -f
	@echo "all needletail images and build cache removed"

## Cargo release build (no Docker)
build:
	cargo build --release -p needletail-server

## Run all core tests
test:
	cargo test -p needletail-core

## Clean cargo build artifacts
clean:
	cargo clean
