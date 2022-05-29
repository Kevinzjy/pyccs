## Developer installation

```
# Clone pyccs & rust-spoa repository
git clone --recursive https://github.com/Kevinzjy/pyccs.git
cd pyccs/misc

# Build docker image
docker build -t manylinux2010_x86_64_rust:v1 ./

# Build wheel
podman run -v \
  /your/path/to/pyccs:/volume:z \
  -w /volume \
  --rm \
  -t localhost/manylinux2010_x86_64_rust:v1 sh ./misc/build-wheels.sh
```