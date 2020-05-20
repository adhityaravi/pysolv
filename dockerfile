#! /bin/sh

# error handling
set -e -x

# running the docker file needs root access. if not root raise error
if [ "$EUID" -ne 0 ]
then
  echo "Root access is necessary. Please run as root"
  exit
fi

# build platform information
# docker image name and tag
DOCKER_IMAGE_NAME="quay.io/pypa/manylinux1_x86_64"
DOCKER_IMAGE_TAG="latest"
# container name
CONTAINER_NAME="pysolv_container"


# function to pull a docker image with the mentioned docker image name and tag
pull_docker() {
  docker pull "$DOCKER_IMAGE_NAME":"$DOCKER_IMAGE_TAG"
}

# create a container with the given name from the specified docker image
create_container() {
  # create the new container
  docker run --name "$CONTAINER_NAME" -d -it -v $(pwd):/io "$DOCKER_IMAGE_NAME":"$DOCKER_IMAGE_TAG"
  # install the required dependencies for pysolv
  # OpenBLAS
  docker exec -it "$CONTAINER_NAME" yum install openblas-devel
}

# start an already existing container
start_container() {
  docker start "$CONTAINER_NAME"
}

# stop a running container
stop_container() {
  docker stop "$CONTAINER_NAME"
}

# remove an existing container
remove_container() {
  docker rm "$CONTAINER_NAME"
}

# check if the container was previously created with a bind mount volume from the pysolv folder, if not stop the
# container, create the container with proper volume binding
check_io() {
  if ! docker exec -it "$CONTAINER_NAME" test -d /io/pysolv
  then
    stop_container
    remove_container
    create_container
  fi
}

# build the pysolv distribution wheel from inside the container
build_wheel() {
  docker exec -it "$CONTAINER_NAME" bash io/build-wheels.sh
}

# Building the wheel for production through a CentOS 5.11 docker (manylinux)
# check if the docker image already exists in the local machine, if not pull it from the registry, create a new
# container, build the distribution wheel and exit the container
if ! docker image inspect "$DOCKER_IMAGE_NAME":"$DOCKER_IMAGE_TAG" > /dev/null 2>&1
then
  pull_docker
  create_container
  build_wheel
  stop_container

# if the docker image already exists, move on to the container checks
else

  # check if the container was already created
  if docker container inspect "$CONTAINER_NAME" > /dev/null 2>&1
  then

    # check if the container is running, if yes, call the build wheel script through the container and sop the container
    if [ ! "$(docker ps -aq -f status=exited -f name=$CONTAINER_NAME)" ]
    then
      check_io
      build_wheel
      stop_container
    # if the container is not running,  start the container, build the wheel and stop the container
    else
      start_container
      check_io
      build_wheel
      stop_container
    fi

  # if there is no container of the specified name for the specified docker image, create a container, build the wheel,
  # and exit the container
  else
    create_container
    build_wheel
    stop_container
  fi

fi
