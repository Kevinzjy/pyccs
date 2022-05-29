## Developer installation

```
# Clone pyccs & rust-spoa repository
git clone --recursive https://github.com/Kevinzjy/pyccs.git
cd pyccs/misc

# Build docker image
docker build -t manylinux2010_x86_64_rust:v1 ./

# Tag image
podman tag 0498fbf5d048 kevinzjy/manylinux2010_x86_64_rust:v1

# Upload image
podman push kevinzjy/manylinux2010_x86_64_rust:v1

# Build wheel
podman run -v \
  /your/path/to/pyccs:/volume:z \
  -w /volume \
  --rm \
  -t kevinzjy/manylinux2010_x86_64_rust:v1 sh ./misc/build-wheels.sh

# Commit
git add --all .
git commit -m "message"
git tag -d v1.1.0 && git tag v1.1.0
git push --tags --force
```